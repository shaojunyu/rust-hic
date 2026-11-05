# Building and development

## Portable Linux build (x86-64)

The supported release target is x86-64 Linux. The release script builds one
statically linked executable, including BWA-MEM3. The build machine needs a
C/C++ compiler and NASM; users running the binary do not. Built-in BWA-MEM3
requires SSE4.2; external Chromap, BWA-MEM and BWA-MEM2 are separate programs
with their own requirements. CPU, kernel and filesystem compatibility still
need testing on the intended machines.

For development, `cargo build --locked` uses the system toolchain and requires
compatible native dependencies. HDF5 is built statically, while other
system-library linkage may depend on the environment; use the release build
below when a portable binary is needed. The checked-in Rust toolchain selects Rust
1.90.0, also the declared minimum supported Rust version; older versions have
not been tested. BWA-MEM3 is always built in; no Cargo feature is needed to
enable it.

For the portable release, run `scripts/build-linux.sh /path/to/new-output` on
an x86-64 system with glibc 2.28, Rust, GCC/G++, CMake, make, NASM, Git, binutils,
and static glibc, libgcc, libstdc++ and zlib development libraries. The script
checks its prerequisites, uses an isolated persistent cache at
`target/release-build-cache` for downloads and build artifacts, patches the
pinned rapidgzip-sys build settings, and rejects binaries with dynamic library
dependencies. The release profile uses Thin LTO, one codegen unit and removes
symbols; this may make debugging harder and increases build time. It copies
license notices alongside the executable. Compiler
requirements for a direct build depend on native dependencies. BWA-MEM3
recommends GCC 15 or Clang 19 for representative performance; older versions
may work but are not guaranteed.

To use a different cache directory, set `RUST_HIC_BUILD_CACHE`. Delete that
directory to force a clean rebuild. A changed Rust/C++ compiler, CMake, NASM,
glibc version, or rapidgzip patch uses a separate cache automatically.

Docker is optional. To exercise the release build from a newer Linux host using
the Rocky Linux 8 / GCC 15 environment defined in this repository:

```sh
docker build -f scripts/Dockerfile.build -t rust-hic-build .
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/src" -w /src rust-hic-build \
  scripts/build-linux.sh /src/release-linux
```

The output directory must not already exist. The image and OS packages are not yet pinned to immutable
versions, so this provides a repeatable build procedure, not a guarantee of
bit-for-bit reproducible binaries. Run `scripts/check.sh` in the build
environment for formatting, Clippy (warnings are errors), and Rust unit tests.
Use `cargo fmt --all` to apply formatting. The Python core tests are described
below.

## GitHub release workflow

`.github/workflows/release.yml` builds the statically linked x86-64 Linux
executable in the repository's Rocky Linux 8 / GCC 15 Docker environment. A
manual workflow run builds and uploads an Actions artifact without publishing.
Pushing a version tag such as `v0.1.0` builds, runs Rust unit tests and
smoke-tests the binary, then creates a GitHub Release with one archive containing
the executable, licenses and documentation. The tag must match the version in
`Cargo.toml`; tags with a prerelease suffix create prereleases. The workflow
needs Actions enabled and repository permission to write release contents.

The workflow caches Cargo downloads and build artifacts. Its Docker base image
and OS packages are not pinned to immutable digests, so binaries are not
guaranteed to be bit-for-bit reproducible. A hosted Actions run is still needed
to validate GitHub-specific permissions and publication.

## Tests

The core test suite uses a small, reproducible dataset in `tests/data` and
covers paired FASTQ indexing, built-in bwa-mem3 mapping against BWA-MEM2, and
pairs cleaning followed by multi-resolution matrix creation. It needs Python 3.9
or newer, `h5py`, and a `bwa-mem2` executable. Build rust-hic first, then run:

```sh
python3 tests/test_core.py
```

If bwa-mem2 is not on `PATH`, set `RUST_HIC_BWA_MEM2` to its executable.
Set `RUST_HIC_BIN` to test a different rust-hic binary. The dataset can
be regenerated with `python3 tests/data/generate.py`; see `tests/README.md`.

## Source layout

- `src/main.rs`: executable entry point and command dispatch.
- `src/cli`: command-line interfaces for index, map, clean and matrix.
- `src/fastq_index`: manifest, index construction, command dispatch, and paired
  chunk streaming in separate modules.
- `src/mapping`: aligner configuration, external process orchestration,
  built-in alignment, and alignment-to-pairs conversion.
- `src/workflow`: YAML workflow configuration and stage orchestration.
- `src/pairs`: pairs sorting, merging, filtering, deduplication and statistics.
- `src/matrix`: binning pairs and creating multi-resolution matrices.

The processing stages run as internal invocations of the same executable; they
are not separate programs to install.
