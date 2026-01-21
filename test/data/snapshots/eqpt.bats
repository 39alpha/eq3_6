setup_file() {
  load '../../test_helper/common-setup'

  _create_snapshot_dir eqpt
}

setup() {
  load '../../test_helper/common-setup'

  _snapshot_setup
}

snapshot() {
  cp "${BATS_TEST_DIRNAME}/../eqpt/${1}.d0" .

  ./eqpt "${1}.d0"

  for ext in d0 d1 d1f s po; do
    cp "${1}.${ext}" "${SNAPSHOT_DIR}/${1}.${ext}"
  done
}

@test "Snapshot (1kb)" {
  snapshot 1kb
}
@test "Snapshot (2kb)" {
  snapshot 2kb
}
@test "Snapshot (500)" {
  snapshot 500
}
@test "Snapshot (5kb)" {
  snapshot 5kb
}
@test "Snapshot (cmp)" {
  snapshot cmp
}
@test "Snapshot (fmt)" {
  snapshot fmt
}
@test "Snapshot (geo)" {
  snapshot geo
}
@test "Snapshot (hmw)" {
  snapshot hmw
}
@test "Snapshot (nh4)" {
  snapshot nh4
}
@test "Snapshot (ph5)" {
  snapshot ph5
}
@test "Snapshot (shv)" {
  snapshot shv
}
@test "Snapshot (sup)" {
  snapshot sup
}
@test "Snapshot (ymp)" {
  snapshot ymp
}
@test "Snapshot (ypf)" {
  snapshot ypf
}
