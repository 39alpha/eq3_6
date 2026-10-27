setup_file() {
  load 'test_helper/common-setup'
  _suite_setup
}

setup() {
  load 'test_helper/common-setup'
  _common_setup
}

function is_setup {
  run which eqpt
  assert_output "${BATS_TEST_TMPDIR}/bin/eqpt"

  run pwd
  assert_output "${BATS_TEST_TMPDIR}"
}

@test "Binary Exists" {
  is_setup
}

@test "Correctly replaces extension" {
  is_setup

  cp "${EQ36_TEST_DATA_DIR}/eqpt/1kb.d0" 1kb.d0
  eqpt 1kb.d0
  assert_exists 1kb.d1
  assert_exists 1kb.d1f
  assert_exists 1kb.po
  assert_exists 1kb.s
}

@test "Actual extension is not relevant" {
  is_setup

  cp "${EQ36_TEST_DATA_DIR}/eqpt/1kb.d0" data0.1kb
  eqpt data0.1kb
  assert_exists data0.d1
  assert_exists data0.d1f
  assert_exists data0.po
  assert_exists data0.s
}

@test "Handles no extension" {
  is_setup

  cp "${EQ36_TEST_DATA_DIR}/eqpt/1kb.d0" 1kb
  eqpt 1kb
  assert_exists 1kb.d1
  assert_exists 1kb.d1f
  assert_exists 1kb.po
  assert_exists 1kb.s
}

check_output() {
  is_setup

  local -r PROBLEM="${1}"

  cp "${EQ36_TEST_DATA_DIR}/eqpt/${PROBLEM}.d0" .
  for ext in d1 d1f s po; do
    cp "${EQ36_TEST_DATA_DIR}/eqpt/${PROBLEM}.${ext}" "expected.${ext}"
  done

  eqpt "${PROBLEM}.d0"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.po
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.po

  for ext in d1f s po; do
    assert_files_almost_same "${SNAPSHOT_DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}"
    assert_files_equal "${PROBLEM}.${ext}" "expected.${ext}"
  done
  assert_files_equal "${PROBLEM}.d1" "expected.d1"
}

declare -a problems=(1kb 2kb 500 5kb cmp fmt geo hmw nh4 ph5 shv sup ymp ypf)

for problem in "${problems[@]}"; do
  bats_test_function --description "Correct output (${problem})" -- check_output "${problem}"
done
