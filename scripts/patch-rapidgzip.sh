#!/usr/bin/env bash
# Patch the fetched rapidgzip-sys 1.2.1 source in the active Cargo cache.
set -euo pipefail
native=("${CARGO_HOME:-$HOME/.cargo}"/registry/src/*/rapidgzip-sys-1.2.1)
[[ ${#native[@]} == 1 && -d ${native[0]} ]] || { echo 'rapidgzip-sys source not found' >&2; exit 1; }
cmake=${native[0]}/native/rapidgzip-capi/CMakeLists.txt
archive=${native[0]}/vendor/librapidarchive/CMakeLists.txt
if ! grep -q 'set(LIBRAPIDARCHIVE_WITH_RPMALLOC OFF CACHE BOOL "" FORCE)' "$cmake"; then
    sed -i 's/set(LIBRAPIDARCHIVE_WITH_RPMALLOC ON CACHE BOOL "" FORCE)/set(LIBRAPIDARCHIVE_WITH_RPMALLOC OFF CACHE BOOL "" FORCE)/' "$cmake"
fi
if ! grep -q 'set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)' "$archive"; then
    sed -i 's/set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)/set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)/' "$archive"
fi
grep -q 'set(LIBRAPIDARCHIVE_WITH_RPMALLOC OFF CACHE BOOL "" FORCE)' "$cmake"
grep -q 'set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)' "$archive"
# The native wrapper otherwise explicitly requests the shared C++ runtime.
if ! grep -q 'static=stdc++' "${native[0]}/build.rs"; then
    sed -i 's/dylib=stdc++/static=stdc++/' "${native[0]}/build.rs"
fi
grep -q 'static=stdc++' "${native[0]}/build.rs"
