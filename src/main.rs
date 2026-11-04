//! Command dispatch for the single rust-hic executable.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
mod cli;
mod fastq_index;
mod mapping;
mod matrix;
mod pairs;
mod stage;
mod workflow;

use clap::{CommandFactory, Parser};
use std::io::Write;

fn main() {
    let args: Vec<_> = std::env::args_os().collect();
    if args.get(1).is_some_and(|arg| arg == "index")
        && args.get(2).is_some_and(|arg| arg == "--internal")
    {
        if let Err(e) = fastq_index::main() {
            eprintln!("error: {e:#}");
            std::process::exit(1);
        }
        return;
    }
    let invoked = args.first().cloned().unwrap_or_default();
    if workflow::stage::dispatch(&invoked) {
        return;
    }
    let cli = cli::Cli::parse();
    let result = match cli.command {
        Some(cli::Command::ExampleConfig) => std::io::stdout()
            .lock()
            .write_all(include_bytes!("../docs/example.config.yaml"))
            .map_err(Into::into),
        Some(cli::Command::Run { config, check }) => workflow::run(&config, check),
        Some(cli::Command::Index(args)) => cli::index::run(args),
        Some(cli::Command::Map(args)) => cli::map::run(args),
        Some(cli::Command::Clean(args)) => cli::postprocess::clean(args),
        Some(cli::Command::Matrix(args)) => cli::postprocess::matrix(args),
        Some(cli::Command::Version) => {
            println!("rust-hic {}", env!("CARGO_PKG_VERSION"));
            return;
        }
        None => {
            let _ = cli::Cli::command().print_help();
            println!();
            return;
        }
    };
    if let Err(e) = result {
        eprintln!("error: {e:#}");
        std::process::exit(1)
    }
}
