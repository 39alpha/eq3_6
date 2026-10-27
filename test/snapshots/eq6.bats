setup_file() {
  load '../test_helper/common-setup'
  _suite_setup release
}

setup() {
  load '../test_helper/common-setup'
  _common_setup release
}

function is_setup {
  assert_not_equal "${EQ36_TEST_SNAPSHOTS_DIR}" ""

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

check_output() {
  is_setup

  local -r DIR="${1}"
  local -r PROBLEM="${2}"
  local -r SNAPSHOT_DIR="${EQ36_TEST_SNAPSHOTS_DIR}/eq6"
  shift 2

  cp "${EQ36_TEST_DATA_DIR}/eqpt/${DIR}.d0" .
  cp "${EQ36_TEST_DATA_DIR}/eq6/${DIR}/${PROBLEM}.6i" .

  run eqpt "${DIR}.d0"
  run eq6 "${DIR}.d1" "${PROBLEM}.6i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.6o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.6o

  for ext in 6o 6p 6tx; do
    assert_files_almost_same "${SNAPSHOT_DIR}/${DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}" 0 0
    assert_files_equal "${SNAPSHOT_DIR}/${DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}"
  done
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
