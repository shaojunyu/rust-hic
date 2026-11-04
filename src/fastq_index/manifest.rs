//! Index manifest and input fingerprints.
use anyhow::Result;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs::{self, File},
    io::{Read, Seek, SeekFrom},
    path::Path,
    time::UNIX_EPOCH,
};

#[derive(Serialize, Deserialize)]
pub(super) struct Input {
    pub(super) path: String,
    pub(super) size: u64,
    pub(super) modified_ns: u128,
    pub(super) edge_sha256: String,
    pub(super) index: String,
    pub(super) uncompressed_size: u64,
}
#[derive(Clone, Serialize, Deserialize)]
pub(super) struct Checkpoint {
    pub(super) compressed_byte_offset: u64,
    pub(super) compressed_bits_to_prime: u8,
    pub(super) uncompressed_offset: u64,
    pub(super) has_window: bool,
}
#[derive(Clone, Serialize, Deserialize)]
pub(super) struct Boundary {
    pub(super) pair: Option<u64>,
    pub(super) r1_uncompressed_offset: u64,
    pub(super) r2_uncompressed_offset: u64,
    pub(super) read_id: String,
    pub(super) r1_checkpoint: Option<Checkpoint>,
    pub(super) r2_checkpoint: Option<Checkpoint>,
}
#[derive(Serialize, Deserialize)]
pub(super) struct Manifest {
    pub(super) version: u32,
    pub(super) description: String,
    pub(super) inputs: Vec<Input>,
    pub(super) pairs: Option<u64>,
    pub(super) chunks: usize,
    pub(super) boundaries: Vec<Boundary>,
}

pub(super) fn fingerprint(path: &Path, index: String, uncompressed_size: u64) -> Result<Input> {
    let path = fs::canonicalize(path)?;
    let mut file = File::open(&path)?;
    let meta = file.metadata()?;
    let size = meta.len();
    let mut h = Sha256::new();
    let mut buf = vec![0; 65536];
    let n = file.read(&mut buf)?;
    h.update(&buf[..n]);
    if size > 65536 {
        file.seek(SeekFrom::Start(size - 65536))?;
        let n = file.read(&mut buf)?;
        h.update(&buf[..n]);
    }
    Ok(Input {
        path: path.to_string_lossy().into_owned(),
        size,
        modified_ns: meta.modified()?.duration_since(UNIX_EPOCH)?.as_nanos(),
        edge_sha256: format!("{:x}", h.finalize()),
        index,
        uncompressed_size,
    })
}
