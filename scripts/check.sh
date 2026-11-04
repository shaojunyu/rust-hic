#!/usr/bin/env bash
# Run local Rust formatting, lint and unit checks with the configured toolchain.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
cd "$root"
# The release cache patches rapidgzip-sys to link libstdc++ statically. GCC
# toolsets can place that archive outside rustc's default library search path.
cxx_archive=$("${CXX:-g++}" -print-file-name=libstdc++.a)
if [[ -f "$cxx_archive" ]]; then
    export RUSTFLAGS="${RUSTFLAGS:+$RUSTFLAGS }-Lnative=$(dirname "$cxx_archive")"
fi
cargo fmt --all -- --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked
