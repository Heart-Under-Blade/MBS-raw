#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=${TMPDIR:-/tmp}/mbs_raw_facet_assignment_test
mkdir -p "$build_dir"

${CXX:-g++} -std=gnu++11 -msse4.2 \
  -I"$repo/src" \
  -I"$repo/src/particle" \
  -I"$repo/src/math" \
  -I"$repo/src/geometry" \
  "$repo/tests/test_facet_assignment.cpp" \
  "$repo/src/geometry/Facet.cpp" \
  "$repo/src/geometry/Polygon.cpp" \
  "$repo/src/geometry/geometry_lib.cpp" \
  -o "$build_dir/test_facet_assignment"

"$build_dir/test_facet_assignment"
