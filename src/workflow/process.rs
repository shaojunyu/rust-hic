//! Launch and supervise internal pipeline stages.
use anyhow::{bail, Context, Result};
use std::{
    fs::{self, File, OpenOptions},
    io,
    os::unix::process::CommandExt,
    path::Path,
    process::{Child, Command, Stdio},
};

pub(crate) fn current_affinity() -> Result<libc::cpu_set_t> {
    let mut mask: libc::cpu_set_t = unsafe { std::mem::zeroed() };
    if unsafe { libc::sched_getaffinity(0, std::mem::size_of_val(&mask), &mut mask) } != 0 {
        return Err(io::Error::last_os_error()).context("read CPU affinity");
    }
    Ok(mask)
}
pub(crate) fn affinity_mask(cpus: Option<&[usize]>) -> Result<Option<libc::cpu_set_t>> {
    let Some(cpus) = cpus else { return Ok(None) };
    if cpus.is_empty() {
        bail!("postprocess_cpus cannot be empty")
    }
    let available = current_affinity()?;
    let mut mask: libc::cpu_set_t = unsafe { std::mem::zeroed() };
    let mut seen = std::collections::HashSet::new();
    for &cpu in cpus {
        if cpu >= libc::CPU_SETSIZE as usize
            || !unsafe { libc::CPU_ISSET(cpu, &available) }
            || !seen.insert(cpu)
        {
            bail!("postprocess CPU {cpu} is duplicated or outside the current CPU affinity")
        }
        unsafe { libc::CPU_SET(cpu, &mut mask) };
    }
    Ok(Some(mask))
}
pub(crate) fn set_affinity(mask: &libc::cpu_set_t) -> io::Result<()> {
    if unsafe { libc::sched_setaffinity(0, std::mem::size_of_val(mask), mask) } == 0 {
        Ok(())
    } else {
        Err(io::Error::last_os_error())
    }
}
pub(super) fn limit_command(cmd: &mut Command, mask: Option<&libc::cpu_set_t>) {
    if let Some(mask) = mask {
        let mask = *mask;
        // sched_setaffinity is async-signal-safe and children inherit it.
        unsafe {
            cmd.pre_exec(move || set_affinity(&mask));
        }
    }
}
pub(crate) fn exec(
    mut cmd: Command,
    dir: &Path,
    name: &str,
    cpus: Option<&libc::cpu_set_t>,
) -> Result<()> {
    let err = File::create(dir.join(format!("{name}.log")))?;
    limit_command(&mut cmd, cpus);
    let status = cmd
        .stdout(Stdio::null())
        .stderr(Stdio::from(err))
        .status()
        .with_context(|| format!("start {name}"))?;
    if !status.success() {
        let log_path = dir.join(format!("{name}.log"));
        let log = fs::read_to_string(&log_path).unwrap_or_default();
        let detail = log
            .lines()
            .rev()
            .find(|line| line.starts_with("error:"))
            .unwrap_or("see log for details");
        bail!(
            "{name} failed: {status}: {detail}; inspect {}",
            log_path.display()
        )
    }
    Ok(())
}
fn stop(children: &mut [Child]) {
    for child in children.iter_mut() {
        if child.try_wait().ok().flatten().is_none() {
            // Each stage has its own process group. This also stops nested
            // range-merge workers if a downstream stage fails.
            unsafe {
                libc::kill(-(child.id() as libc::pid_t), libc::SIGKILL);
            }
            let _ = child.kill();
        }
    }
    for child in children.iter_mut() {
        let _ = child.wait();
    }
}
pub(crate) fn pipe(
    mut commands: Vec<(String, Command)>,
    output: &Path,
    logs: &Path,
    cpus: Option<&libc::cpu_set_t>,
) -> Result<()> {
    let mut children: Vec<Child> = Vec::new();
    let mut previous = None;
    let count = commands.len();
    for (idx, (name, cmd)) in commands.iter_mut().enumerate() {
        cmd.process_group(0);
        limit_command(cmd, cpus);
        let err = match File::create(logs.join(format!("{name}.log"))) {
            Ok(v) => v,
            Err(e) => {
                stop(&mut children);
                return Err(e.into());
            }
        };
        cmd.stderr(Stdio::from(err));
        if let Some(input) = previous.take() {
            cmd.stdin(Stdio::from(input));
        } else {
            cmd.stdin(Stdio::null());
        }
        if idx + 1 == count {
            let f = match OpenOptions::new().create_new(true).write(true).open(output) {
                Ok(v) => v,
                Err(e) => {
                    stop(&mut children);
                    return Err(e.into());
                }
            };
            cmd.stdout(Stdio::from(f));
        } else {
            cmd.stdout(Stdio::piped());
        }
        match cmd.spawn() {
            Ok(mut child) => {
                previous = child.stdout.take();
                children.push(child)
            }
            Err(e) => {
                stop(&mut children);
                return Err(e).with_context(|| format!("start {name}"));
            }
        }
    }
    for i in (0..children.len()).rev() {
        match children[i].wait() {
            Ok(status) if status.success() => {}
            Ok(status) => {
                let name = commands[i].0.clone();
                stop(&mut children);
                let details = commands
                    .iter()
                    .filter_map(|(stage, _)| {
                        let content = fs::read_to_string(logs.join(format!("{stage}.log"))).ok()?;
                        let error = content
                            .lines()
                            .rev()
                            .find(|line| line.starts_with("error:"))?;
                        Some(format!("{stage}: {error}"))
                    })
                    .collect::<Vec<_>>()
                    .join("; ");
                bail!(
                    "{name} failed: {status}; {details}; inspect {}",
                    logs.join(format!("{name}.log")).display()
                )
            }
            Err(e) => {
                stop(&mut children);
                return Err(e.into());
            }
        }
    }
    Ok(())
}
