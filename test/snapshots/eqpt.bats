setup_file() {
  load '../test_helper/common-setup'
  _suite_setup
}

setup() {
  load '../test_helper/common-setup'
  _common_setup
}

function is_setup {
  assert_not_equal "${EQ36_TEST_SNAPSHOTS_DIR}" ""

  run which eqpt
  assert_output "${BATS_TEST_TMPDIR}/bin/eqpt"

  run pwd
  assert_output "${BATS_TEST_TMPDIR}"
}

@test "Binary Exists" {
  is_setup
}

check_output() {
  is_setup

  local -r PROBLEM="${1}"
  local -r SNAPSHOT_DIR="${EQ36_TEST_SNAPSHOTS_DIR}/eqpt"

  cp "${EQ36_TEST_DATA_DIR}/eqpt/${PROBLEM}.d0" .

  eqpt "${PROBLEM}.d0"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.po
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.po

  for ext in d1f s po; do
    assert_files_almost_same "${SNAPSHOT_DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}" 0 0
    assert_files_equal "${SNAPSHOT_DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}"
  done
  assert_files_equal "${SNAPSHOT_DIR}/${PROBLEM}.d1" "${PROBLEM}.d1"
}

declare -a problems=(1kb 2kb 500 5kb cmp fmt geo hmw nh4 ph5 shv sup ymp ypf)

for problem in "${problems[@]}"; do
  bats_test_function --description "Correct output (${problem})" -- check_output "${problem}"
done
