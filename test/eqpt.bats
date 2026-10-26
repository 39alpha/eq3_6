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

  cp "${BATS_TEST_DIRNAME}/data/eqpt/1kb.d0" 1kb.d0
  eqpt 1kb.d0
  assert_exists 1kb.d1
  assert_exists 1kb.d1f
  assert_exists 1kb.po
  assert_exists 1kb.s
}

@test "Actual extension is not relevant" {
  is_setup

  cp "${BATS_TEST_DIRNAME}/data/eqpt/1kb.d0" data0.1kb
  eqpt data0.1kb
  assert_exists data0.d1
  assert_exists data0.d1f
  assert_exists data0.po
  assert_exists data0.s
}

@test "Handles no extension" {
  is_setup

  cp "${BATS_TEST_DIRNAME}/data/eqpt/1kb.d0" 1kb
  eqpt 1kb
  assert_exists 1kb.d1
  assert_exists 1kb.d1f
  assert_exists 1kb.po
  assert_exists 1kb.s
}

check_output() {
  is_setup

  local -r PROBLEM="${1}"
  local -r SNAPSHOT_DIR="$(_snapshot_dir eqpt)"

  cp "${BATS_TEST_DIRNAME}/data/eqpt/${PROBLEM}.d0" .
  for ext in d1 d1f s po; do
    cp "${BATS_TEST_DIRNAME}/data/eqpt/${PROBLEM}.${ext}" "expected.${ext}"
  done

  eqpt "${PROBLEM}.d0"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.po
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.po

  for ext in d1 d1f s po; do
    assert_files_equal "${PROBLEM}.${ext}" "expected.${ext}"
  done

  for ext in d1f s po; do
    assert_files_almost_same "${SNAPSHOT_DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}" 0 0
    assert_files_equal "${SNAPSHOT_DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}"
  done
  assert_files_equal "${SNAPSHOT_DIR}/${PROBLEM}.d1" "${PROBLEM}.d1"
}

@test "Correct output (1kb)" {
  check_output 1kb
}
@test "Correct output (2kb)" {
  check_output 2kb
}
@test "Correct output (500)" {
  check_output 500
}
@test "Correct output (5kb)" {
  check_output 5kb
}
@test "Correct output (cmp)" {
  check_output cmp
}
@test "Correct output (fmt)" {
  check_output fmt
}
@test "Correct output (geo)" {
  check_output geo
}
@test "Correct output (hmw)" {
  check_output hmw
}
@test "Correct output (nh4)" {
  check_output nh4
}
@test "Correct output (ph5)" {
  check_output ph5
}
@test "Correct output (shv)" {
  check_output shv
}
@test "Correct output (sup)" {
  check_output sup
}
@test "Correct output (ymp)" {
  check_output ymp
}
@test "Correct output (ypf)" {
  check_output ypf
}
