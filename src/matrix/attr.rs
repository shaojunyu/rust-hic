//! Shared HDF5 attribute helpers.
use anyhow::Result;
use hdf5::{types::VarLenUnicode, Group};
use std::str::FromStr;

pub(super) fn str_attr(group: &Group, key: &str, value: &str) -> Result<()> {
    group
        .new_attr::<VarLenUnicode>()
        .shape(())
        .create(key)?
        .write_scalar(&VarLenUnicode::from_str(value)?)?;
    Ok(())
}

pub(super) fn int_attr(group: &Group, key: &str, value: i64) -> Result<()> {
    group
        .new_attr::<i64>()
        .shape(())
        .create(key)?
        .write_scalar(&value)?;
    Ok(())
}
