mod attr;
pub(crate) mod bin_pairs;
pub(crate) mod multiresolution;

/// Internal setting passed to the matrix balancing stage.
pub(crate) const BALANCE_MAX_ITERS_ENV: &str = "RUST_HIC_BALANCE_MAX_ITERS";
