//! Slurm job-array execution: index, parallel alignment, and final processing.
mod state;
mod status;
mod submit;
mod worker;

pub(crate) use status::{cancel, status};
pub(crate) use submit::submit;
pub(crate) use worker::task;
