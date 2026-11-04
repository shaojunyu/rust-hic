#!/usr/bin/env bash
# Run on an x86-64 RHEL 8 compatible build host with the dependencies in README.
set -euo pipefail
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
out=${1:?usage: build-linux.sh NEW_OUTPUT_DIRECTORY}
if [[ -e "$out" || -L "$out" ]]; then echo "output already exists: $out" >&2; exit 1; fi
out_parent=$(dirname -- "$out")
[[ -d "$out_parent" && -w "$out_parent" ]] || { echo "output parent must be an existing writable directory: $out_parent" >&2; exit 1; }
[[ $(uname -m) == x86_64 ]] || { echo 'x86-64 required' >&2; exit 1; }
[[ $(getconf GNU_LIBC_VERSION) == 'glibc 2.28' ]] || { echo 'build on a glibc 2.28 system (RHEL 8)' >&2; exit 1; }
for tool in cargo rustc "${CC:-gcc}" "${CXX:-g++}" cmake make nasm git ar readelf flock sha256sum; do
    command -v "$tool" >/dev/null || { echo "required build tool not found: $tool" >&2; exit 1; }
done
"${CC:-gcc}" --version | head -n 1
"${CXX:-g++}" --version | head -n 1
rustc --version
cmake --version | head -n 1
nasm -v
cxx_archive=$("${CXX:-g++}" -print-file-name=libstdc++.a)
[[ -f "$cxx_archive" ]] || { echo 'static libstdc++ is required' >&2; exit 1; }
for archive in libc.a libgcc.a libz.a; do
    archive_path=$("${CC:-gcc}" -print-file-name="$archive")
    [[ -f "$archive_path" ]] || { echo "static $archive is required" >&2; exit 1; }
done
cd "$root"
# Avoid inheriting machine-specific HDF5 paths or CPU tuning in release builds.
unset HDF5_DIR RUSTFLAGS CFLAGS CXXFLAGS
cache_root=${RUST_HIC_BUILD_CACHE:-"$root/target/release-build-cache"}
mkdir -p "$cache_root"
cache_root=$(cd -- "$cache_root" && pwd -P)
# Native build results depend on the toolchain and on our rapidgzip patch.
cache_key=$(
    { rustc -vV; cargo -V; "${CC:-gcc}" --version; "${CXX:-g++}" --version;
      command -v "${CC:-gcc}"; command -v "${CXX:-g++}";
      cmake --version; nasm -v; getconf GNU_LIBC_VERSION;
      sha256sum "$root/scripts/patch-rapidgzip.sh"; } | sha256sum | cut -d ' ' -f 1
)
cache_dir="$cache_root/$cache_key"
mkdir -p "$cache_dir"
# The source patch and Cargo target directory must not be updated concurrently.
exec 9>"$cache_dir/build.lock"
flock -x 9
export CARGO_HOME="$cache_dir/cargo" CARGO_TARGET_DIR="$cache_dir/target"
if ! cargo fetch --locked --offline 2>/dev/null; then
    cargo fetch --locked
fi
"$root/scripts/patch-rapidgzip.sh"
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_RUSTFLAGS="-C target-feature=+crt-static -L native=$(dirname "$cxx_archive")"
cargo build --release --locked --target x86_64-unknown-linux-gnu --bin rust-hic
if [[ ${RUST_HIC_TEST_RELEASE:-0} == 1 ]]; then
    cargo test --release --locked --target x86_64-unknown-linux-gnu
fi
binary="$CARGO_TARGET_DIR/x86_64-unknown-linux-gnu/release/rust-hic"
if readelf -l "$binary" | grep -q INTERP || readelf -d "$binary" | grep -q NEEDED; then
    echo 'release binary is dynamically linked' >&2
    exit 1
fi
mkdir -p "$out/licenses"
cp "$binary" "$out/rust-hic"
for pkg in glibc libgcc zlib; do
    if [[ -d /usr/share/licenses/$pkg ]]; then cp -r "/usr/share/licenses/$pkg" "$out/licenses/"; fi
done
cp "$root/licenses/"* "$out/licenses/"
hdf5_src=("$CARGO_HOME"/registry/src/*/hdf5-metno-src-*/ext/hdf5/LICENSE)
if [[ -f ${hdf5_src[0]} ]]; then cp "${hdf5_src[0]}" "$out/licenses/hdf5-LICENSE"; fi
cp "$root/README.md" "$out/README.md"
cp "$root/LICENSE" "$out/LICENSE"
cp -r "$root/docs" "$out/docs"
