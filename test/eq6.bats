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

  run which eq6
  assert_output "${BATS_TEST_TMPDIR}/bin/eq6"

  run pwd
  assert_output "${BATS_TEST_TMPDIR}"
}

@test "Binary Exists" {
  is_setup
}

@test "Correctly replaces extension" {
  is_setup

  cp "${EQ36_TEST_DATA_DIR}/eqpt/cmp.d1" .
  cp "${EQ36_TEST_DATA_DIR}/eq6/cmp/crisqtz.6i" .
  eq6 cmp.d1 crisqtz.6i
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

@test "Actual extension is not relevant" {
  is_setup

  cp "${EQ36_TEST_DATA_DIR}/eqpt/cmp.d1" .
  cp "${EQ36_TEST_DATA_DIR}/eq6/cmp/crisqtz.6i" crisqtz.prob
  eq6 cmp.d1 crisqtz.prob
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

@test "Handles no extension" {
  cp "${EQ36_TEST_DATA_DIR}/eqpt/cmp.d1" .
  cp "${EQ36_TEST_DATA_DIR}/eq6/cmp/crisqtz.6i" crisqtz
  eq6 cmp.d1 crisqtz
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

check_output() {
  is_setup

  local -r DIR="${1}"
  local -r PROBLEM="${2}"
  shift 2

  cp "${EQ36_TEST_DATA_DIR}/eqpt/${DIR}.d0" .
  cp "${EQ36_TEST_DATA_DIR}/eq6/${DIR}/${PROBLEM}.6i" .
  for ext in 6i 6o 6p; do
    cp "${EQ36_TEST_DATA_DIR}/eq6/${DIR}/${PROBLEM}.${ext}" "expected.${ext}"
  done

  run eqpt "${DIR}.d0"
  run eq6 "${DIR}.d1" "${PROBLEM}.6i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.6o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.6o

  assert_files_almost_same  "expected.6p" "${PROBLEM}.6p" 1e-12 3e-5
  assert_files_almost_same  "expected.6tx" "${PROBLEM}.6tx"
}

for thermo in cmp fmt hmw ymp ypf; do
  declare -a problems
  case "${thermo}" in
    "cmp")
      problems=(crisqtz dedolo heatqf heatsw heatswfl j13wsf j13wtitr j13wtuff methane micro microft pptcal pptmins pptqtz pptqtza pyrsw rwssdiag rwtitr swtitr swxrca swxrcaft)
      ;;
    "fmt")
      problems=(c4pgwbN2 f24vc7b3 f24vc7k4 f24vc7m gypnaclx)
      ;;
    "hmw")
      problems=(calhal evapsw evswgyha fwbrmix gypanhy mgso4 swv1sxk)
      ;;
    "ymp")
      problems=(crisqtz dedolo heatqf heatsw j13wsf j13wtitr j13wtuff methane micro microft pptcal pptmins pptqtz pptqtza pyrsw rwssdiag rwtitr swtitr swxrca swxrcaft)
      ;;
    "ypf")
      problems=(calhal90 evapsw60)
      ;;
  esac

  for problem in "${problems[@]}"; do
    bats_test_function --description "Correct output (${thermo}/${problem})" -- check_output "${thermo}" "${problem}"
  done
done
