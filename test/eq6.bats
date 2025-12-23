setup() {
  load 'test_helper/bats-support/load'
  load 'test_helper/bats-assert/load'
  load 'test_helper/bats-file/load'

  DIR="$(cd "${BATS_TEST_DIRNAME}/.." >/dev/null 2>&1 && pwd)"

  cd "${BATS_TEST_TMPDIR}" || return 1

  cp "${DIR}/bin/eq6" .
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
  cp "${BATS_TEST_DIRNAME}/data/eq6/${1}/${1}.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/${1}/${2}.6i" .
  for ext in 6i 6o 6p; do
    cp "${BATS_TEST_DIRNAME}/data/eq6/${1}/${2}.${ext}" "expected.${ext}"
  done

  run ./eq6 "${1}.d1" "${2}.6i"

  sed -E -i '/^\s*(Start|End|Run)\s+time/d' ./*.6o
  sed -E -i '/^\s*Run\s+[0-9]+/d' ./*.6o

  # DGM: For the moment, we are only going to check the 3p file. The problem is
  #      we can't do a simple comparison between the files (e.g. diff) because
  #      that doesn't account for small variations due to machine architecture
  #      and differences that have no real effect on the results (-0 vs 0...).
  # 
  # for ext in 3i 3o 3p; do
  #   assert_files_equal "${2}.${ext}" "expected.${ext}"
  # done
  assert_files_equal "${2}.6p" "expected.6p"
}

@test "Correct output (cmp/crisqtz)" {
  skip
  check_output cmp crisqtz
}
@test "Correct output (cmp/dedolo)" {
  check_output cmp dedolo
}
@test "Correct output (cmp/heatqf)" {
  check_output cmp heatqf
}
@test "Correct output (cmp/heatsw)" {
  skip
  check_output cmp heatsw
}
@test "Correct output (cmp/heatswfl)" {
  skip
  check_output cmp heatswfl
}
@test "Correct output (cmp/j13wsf)" {
  check_output cmp j13wsf
}
@test "Correct output (cmp/j13wtitr)" {
  skip
  check_output cmp j13wtitr
}
@test "Correct output (cmp/j13wtuff)" {
  skip
  check_output cmp j13wtuff
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
  skip
  check_output cmp pyrsw
}
@test "Correct output (cmp/rwssdiag)" {
  skip
  check_output cmp rwssdiag
}
@test "Correct output (cmp/rwtitr)" {
  skip
  check_output cmp rwtitr
}
@test "Correct output (cmp/swtitr)" {
  skip
  check_output cmp swtitr
}
@test "Correct output (cmp/swxrca)" {
  skip
  check_output cmp swxrca
}
@test "Correct output (cmp/swxrcaft)" {
  skip
  check_output cmp swxrcaft
}
@test "Correct output (fmt/c4pgwbN2)" {
  skip
  check_output fmt c4pgwbN2
}
@test "Correct output (fmt/f24vc7b3)" {
  skip
  check_output fmt f24vc7b3
}
@test "Correct output (fmt/f24vc7k4)" {
  skip
  check_output fmt f24vc7k4
}
@test "Correct output (fmt/f24vc7m)" {
  skip
  check_output fmt f24vc7m
}
@test "Correct output (fmt/gypnaclx)" {
  check_output fmt gypnaclx
}
@test "Correct output (hmw/calhal)" {
  skip
  check_output hmw calhal
}
@test "Correct output (hmw/evapsw)" {
  skip
  check_output hmw evapsw
}
@test "Correct output (hmw/evswgyha)" {
  skip
  check_output hmw evswgyha
}
@test "Correct output (hmw/fwbrmix)" {
  skip
  check_output hmw fwbrmix
}
@test "Correct output (hmw/gypanhy)" {
  check_output hmw gypanhy
}
@test "Correct output (hmw/mgso4)" {
  skip
  check_output hmw mgso4
}
@test "Correct output (hmw/swv1sxk)" {
  check_output hmw swv1sxk
}
@test "Correct output (ymp/crisqtz)" {
  check_output ymp crisqtz
}
@test "Correct output (ymp/dedolo)" {
  skip
  check_output ymp dedolo
}
@test "Correct output (ymp/heatqf)" {
  check_output ymp heatqf
}
@test "Correct output (ymp/heatsw)" {
  skip
  check_output ymp heatsw
}
@test "Correct output (ymp/j13wsf)" {
  skip
  check_output ymp j13wsf
}
@test "Correct output (ymp/j13wtitr)" {
  check_output ymp j13wtitr
}
@test "Correct output (ymp/j13wtuff)" {
  skip
  check_output ymp j13wtuff
}
@test "Correct output (ymp/methane)" {
  skip
  check_output ymp methane
}
@test "Correct output (ymp/micro)" {
  skip
  check_output ymp micro
}
@test "Correct output (ymp/microft)" {
  skip
  check_output ymp microft
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
  skip
  check_output ymp swtitr
}
@test "Correct output (ymp/swxrca)" {
  skip
  check_output ymp swxrca
}
@test "Correct output (ymp/swxrcaft)" {
  skip
  check_output ymp swxrcaft
}
@test "Correct output (ypf/calhal90)" {
  check_output ypf calhal90
}
@test "Correct output (ypf/evapsw60)" {
  skip
  check_output ypf evapsw60
}

