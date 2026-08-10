#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
build_dir=$(mktemp -d)
trap 'rm -rf "$build_dir"' EXIT

${CXX:-g++} -std=gnu++11 -msse4.2 \
  -I"$repo/src" \
  -I"$repo/src/particle" \
  -I"$repo/src/math" \
  -I"$repo/src/geometry" \
  -I"$repo/src/common" \
  "$repo/tests/test_particle_visibility_reset.cpp" \
  "$repo/src/particle/Particle.cpp" \
  "$repo/src/particle/Hexagonal.cpp" \
  "$repo/src/particle/ConcaveHexagonal.cpp" \
  "$repo/src/particle/Bullet.cpp" \
  "$repo/src/particle/HexagonalAggregate.cpp" \
  "$repo/src/geometry/Facet.cpp" \
  "$repo/src/geometry/Polygon.cpp" \
  "$repo/src/geometry/geometry_lib.cpp" \
  "$repo/src/math/compl.cpp" \
  "$repo/src/common/common.cpp" \
  -o "$build_dir/test_particle_visibility_reset"

"$build_dir/test_particle_visibility_reset"
