setup() {
  load 'test_helper/common-setup'
  _common_setup

  cp "${TEST_TMPDIR}/bin/eqpt" .
}

@test "Binary Exists" {
  assert_exists eqpt
}

@test "Correctly replaces extension" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/1kb.d0" 1kb.d0
  ./eqpt 1kb.d0
  assert_exists 1kb.d1
  assert_exists 1kb.d1f
  assert_exists 1kb.po
  assert_exists 1kb.s
}

@test "Actual extension is not relevant" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/1kb.d0" data0.1kb
  ./eqpt data0.1kb
  assert_exists data0.d1
  assert_exists data0.d1f
  assert_exists data0.po
  assert_exists data0.s
}

@test "Handles no extension" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/1kb.d0" 1kb
  ./eqpt 1kb
  assert_exists 1kb.d1
  assert_exists 1kb.d1f
  assert_exists 1kb.po
  assert_exists 1kb.s
}

check_output() {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/${1}.d0" .
  for ext in d1 d1f s po; do
    cp "${BATS_TEST_DIRNAME}/data/eqpt/${1}.${ext}" "expected.${ext}"
  done

  run ./eqpt "${1}.d0"

  sed -E -i '/^\s*(Start|End|Run)\s+time/d' ./*.po
  sed -E -i '/^\s*Run\s+[0-9]+/d' ./*.po

  for ext in d1 d1f s po; do
    assert_files_equal "${1}.${ext}" "expected.${ext}"
  done
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
