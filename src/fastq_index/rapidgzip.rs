//! Rapidgzip reader for fast indexing. Owns its handle and streams the index to
//! the sparse writer without depending on `rapidgzip::Reader`'s private layout.
use rapidgzip_sys as ffi;
use std::{
    ffi::{c_void, CStr, CString},
    io::{self, Read, Seek, SeekFrom, Write},
    os::unix::ffi::OsStrExt,
    path::Path,
    ptr::NonNull,
};

pub struct RapidgzipReader(NonNull<ffi::rgz_reader_t>);
// The handle is uniquely owned and is never used concurrently.
unsafe impl Send for RapidgzipReader {}

fn message(ptr: *const std::ffi::c_char) -> String {
    if ptr.is_null() {
        "rapidgzip operation failed".into()
    } else {
        // SAFETY: rapidgzip owns this NUL-terminated message until its next call.
        unsafe { CStr::from_ptr(ptr) }
            .to_string_lossy()
            .into_owned()
    }
}
impl RapidgzipReader {
    pub fn open(path: &Path, threads: u32) -> io::Result<Self> {
        let path = CString::new(path.as_os_str().as_bytes())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "path contains NUL"))?;
        let mut config = ffi::rgz_config_t {
            struct_size: 0,
            flags: 0,
            parallelism: 0,
            reserved0: 0,
            chunk_size: 0,
            reserved1: 0,
        };
        // SAFETY: config is valid and path remains alive through the call.
        unsafe { ffi::rgz_config_init(&mut config) };
        config.parallelism = threads;
        config.chunk_size = 4 * 1024 * 1024;
        let mut status = ffi::rgz_status_t::RGZ_STATUS_OK;
        let ptr = unsafe { ffi::rgz_open_path_ex(path.as_ptr(), &config, &mut status) };
        let ptr = NonNull::new(ptr).ok_or_else(|| {
            io::Error::other(format!(
                "rapidgzip open ({status:?}): {}",
                message(unsafe { ffi::rgz_last_global_error() })
            ))
        })?;
        let mut reader = Self(ptr);
        let status = unsafe { ffi::rgz_set_keep_index(reader.0.as_ptr(), true) };
        reader.check(status)?;
        Ok(reader)
    }
    pub fn import_index(&mut self, path: &Path) -> io::Result<()> {
        let path = CString::new(path.as_os_str().as_bytes())
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "index path contains NUL"))?;
        // SAFETY: the owned handle is valid and path remains live during the call.
        let status = unsafe { ffi::rgz_import_index_path(self.0.as_ptr(), path.as_ptr()) };
        self.check(status)
    }
    fn check(&mut self, status: ffi::rgz_status_t) -> io::Result<()> {
        if status == ffi::rgz_status_t::RGZ_STATUS_OK {
            Ok(())
        } else {
            Err(io::Error::other(format!(
                "rapidgzip {status:?}: {}",
                message(unsafe { ffi::rgz_last_error(self.0.as_ptr()) })
            )))
        }
    }
    pub fn export_to(&mut self, writer: &mut dyn Write) -> io::Result<()> {
        struct State<'a> {
            writer: &'a mut dyn Write,
            error: Option<io::Error>,
            panicked: bool,
        }
        unsafe extern "C" fn write_cb(
            user: *mut c_void,
            buffer: *const c_void,
            size: usize,
        ) -> usize {
            // SAFETY: user points to the live State for the synchronous export call.
            let state = &mut *(user as *mut State<'_>);
            if size == 0 {
                return 0;
            }
            if buffer.is_null() {
                state.error = Some(io::Error::other("null rapidgzip output buffer"));
                return 0;
            }
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let bytes = std::slice::from_raw_parts(buffer.cast::<u8>(), size);
                state.writer.write_all(bytes)
            }));
            match result {
                Ok(Ok(())) => size,
                Ok(Err(err)) => {
                    state.error = Some(err);
                    0
                }
                Err(_) => {
                    state.panicked = true;
                    0
                }
            }
        }
        let mut state = State {
            writer,
            error: None,
            panicked: false,
        };
        // SAFETY: the native callback is synchronous and state stays alive throughout.
        let status = unsafe {
            ffi::rgz_export_index_write(
                self.0.as_ptr(),
                write_cb,
                (&mut state as *mut State<'_>).cast(),
            )
        };
        if let Some(err) = state.error {
            return Err(err);
        }
        if state.panicked {
            return Err(io::Error::other("index writer panicked"));
        }
        self.check(status)
    }
}
impl Drop for RapidgzipReader {
    fn drop(&mut self) {
        unsafe { ffi::rgz_close(self.0.as_ptr()) }
    }
}
impl Read for RapidgzipReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        let mut count = 0;
        let status = unsafe {
            ffi::rgz_read(
                self.0.as_ptr(),
                buf.as_mut_ptr().cast(),
                buf.len(),
                &mut count,
            )
        };
        if status != ffi::rgz_status_t::RGZ_STATUS_EOF {
            self.check(status)?;
        }
        Ok(count)
    }
}
impl Seek for RapidgzipReader {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let offset = match pos {
            SeekFrom::Start(n) => n,
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "only absolute seeks are supported",
                ))
            }
        };
        let status = unsafe { ffi::rgz_seek_to(self.0.as_ptr(), offset) };
        self.check(status)?;
        Ok(offset)
    }
}
