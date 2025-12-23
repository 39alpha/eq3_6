setup() {
  load 'test_helper/common-setup'
  _common_setup

  cp "${TEST_TMPDIR}/bin/eq6" .
}

@test "Binary Exists" {
  assert_exists eq6
}

@test "Correctly replaces extension" {
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/crisqtz.6i" .
  ./eq6 cmp.d1 crisqtz.6i
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

@test "Actual extension is not relevant" {
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/crisqtz.6i" crisqtz.prob
  ./eq6 cmp.d1 crisqtz.prob
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

@test "Handles no extension" {
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/crisqtz.6i" crisqtz
  ./eq6 cmp.d1 crisqtz
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

check_output() {
  local -r DIR="${1}"
  local -r PROBLEM="${2}"
  shift 2

  cp "${BATS_TEST_DIRNAME}/data/eq6/${DIR}/${DIR}.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/${DIR}/${PROBLEM}.6i" .
  for ext in 6i 6o 6p; do
    cp "${BATS_TEST_DIRNAME}/data/eq6/${DIR}/${PROBLEM}.${ext}" "expected.${ext}"
  done

  run ./eq6 "${DIR}.d1" "${PROBLEM}.6i"

  # DGM: We are not going to check the 6o files for a couple of reasons:
  #        1. They can differ in the number of iterations taken. While that
  #           matters, it also kinda doesn't...
  #        2. We are checking the pickup, 6tx and csv, which have all of the
  #           the majority of the computed values of note.
  #
  # sed -E -i '/^\s*(Start|End|Run)\s+time/d' ./*.6o
  # sed -E -i '/^\s*Run\s+[0-9]+/d' ./*.6o

  for ext in 6p 6tx csv; do
    assert_files_almost_same  "expected.${ext}" "${PROBLEM}.${ext}" "${@}"
  done
}

@test "Correct output (cmp/crisqtz)" {
  check_output cmp crisqtz 1e-14
}
@test "Correct output (cmp/dedolo)" {
  check_output cmp dedolo
}
@test "Correct output (cmp/heatqf)" {
  check_output cmp heatqf
}
@test "Correct output (cmp/heatsw)" {
  check_output cmp heatsw 1e-15
}
@test "Correct output (cmp/heatswfl)" {
  check_output cmp heatswfl 1e-6
}
@test "Correct output (cmp/j13wsf)" {
  check_output cmp j13wsf
}
@test "Correct output (cmp/j13wtitr)" {
  check_output cmp j13wtitr 1e-14
}
@test "Correct output (cmp/j13wtuff)" {
  check_output cmp j13wtuff 1e-9
}
@test "Correct output (cmp/methane)" {
  check_output cmp methane
}
@test "Correct output (cmp/micro)" {
  check_output cmp micro
}
@test "Correct output (cmp/microft)" {
  check_output cmp microft
}
@test "Correct output (cmp/pptcal)" {
  check_output cmp pptcal
}
@test "Correct output (cmp/pptmins)" {
  check_output cmp pptmins
}
@test "Correct output (cmp/pptqtz)" {
  check_output cmp pptqtz
}
@test "Correct output (cmp/pptqtza)" {
  check_output cmp pptqtza
}
@test "Correct output (cmp/pyrsw)" {
  check_output cmp pyrsw 1e-19
}
@test "Correct output (cmp/rwssdiag)" {
  check_output cmp rwssdiag 1e-8
}
@test "Correct output (cmp/rwtitr)" {
  check_output cmp rwtitr 1e-17
}
@test "Correct output (cmp/swtitr)" {
  check_output cmp swtitr 1e-9
}
@test "Correct output (cmp/swxrca)" {
  check_output cmp swxrca 1e-6
}
@test "Correct output (cmp/swxrcaft)" {
  check_output cmp swxrcaft 1e-5
}
@test "Correct output (fmt/c4pgwbN2)" {
  check_output fmt c4pgwbN2 1e-11
}
@test "Correct output (fmt/f24vc7b3)" {
  check_output fmt f24vc7b3 1e-13
}
@test "Correct output (fmt/f24vc7k4)" {
  check_output fmt f24vc7k4 1e-11
}
@test "Correct output (fmt/f24vc7m)" {
  check_output fmt f24vc7m 1e-13
}
@test "Correct output (fmt/gypnaclx)" {
  check_output fmt gypnaclx
}
@test "Correct output (hmw/calhal)" {
  check_output hmw calhal 1e-14
}
@test "Correct output (hmw/evapsw)" {
  check_output hmw evapsw 1e-7
}
@test "Correct output (hmw/evswgyha)" {
  check_output hmw evswgyha 1e-12
}
@test "Correct output (hmw/fwbrmix)" {
  check_output hmw fwbrmix 1e-14
}
@test "Correct output (hmw/gypanhy)" {
  check_output hmw gypanhy
}
@test "Correct output (hmw/mgso4)" {
  check_output hmw mgso4 1e-8
}
@test "Correct output (hmw/swv1sxk)" {
  check_output hmw swv1sxk
}
@test "Correct output (ymp/crisqtz)" {
  check_output ymp crisqtz
}
@test "Correct output (ymp/dedolo)" {
  check_output ymp dedolo 1e-13
}
@test "Correct output (ymp/heatqf)" {
  check_output ymp heatqf
}
@test "Correct output (ymp/heatsw)" {
  check_output ymp heatsw 1e-16
}
@test "Correct output (ymp/j13wsf)" {
  check_output ymp j13wsf 1e-9
}
@test "Correct output (ymp/j13wtitr)" {
  check_output ymp j13wtitr
}
@test "Correct output (ymp/j13wtuff)" {
  check_output ymp j13wtuff 1e-6
}
@test "Correct output (ymp/methane)" {
  check_output ymp methane 1e-5
}
@test "Correct output (ymp/micro)" {
  check_output ymp micro 1e-13
}
@test "Correct output (ymp/microft)" {
  check_output ymp microft 1e-4
}
@test "Correct output (ymp/pptcal)" {
  check_output ymp pptcal
}
@test "Correct output (ymp/pptmins)" {
  check_output ymp pptmins
}
@test "Correct output (ymp/pptqtz)" {
  check_output ymp pptqtz
}
@test "Correct output (ymp/pptqtza)" {
  check_output ymp pptqtza
}
@test "Correct output (ymp/pyrsw)" {
  check_output ymp pyrsw
}
@test "Correct output (ymp/rwssdiag)" {
  check_output ymp rwssdiag
}
@test "Correct output (ymp/rwtitr)" {
  check_output ymp rwtitr
}
@test "Correct output (ymp/swtitr)" {
  check_output ymp swtitr 1e-13
}
@test "Correct output (ymp/swxrca)" {
  check_output ymp swxrca 1e-13
}
@test "Correct output (ymp/swxrcaft)" {
  check_output ymp swxrcaft 1e-8
}
@test "Correct output (ypf/calhal90)" {
  check_output ypf calhal90
}
@test "Correct output (ypf/evapsw60)" {
  check_output ypf evapsw60 1e-12
}
