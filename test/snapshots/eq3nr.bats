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

  run which eq3nr
  assert_output "${BATS_TEST_TMPDIR}/bin/eq3nr"

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
  local -r SNAPSHOT_DIR="${EQ36_TEST_SNAPSHOTS_DIR}/eq3nr"
  shift 2

  cp "${EQ36_TEST_DATA_DIR}/eqpt/${DIR}.d0" .
  cp "${EQ36_TEST_DATA_DIR}/eq3nr/${DIR}/${PROBLEM}.3i" .

  run eqpt "${DIR}.d0"
  run eq3nr "${DIR}.d1" "${PROBLEM}.3i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.3o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.3o

  for ext in 3o 3p; do
    assert_files_almost_same "${SNAPSHOT_DIR}/${DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}" 0 0
    assert_files_equal "${SNAPSHOT_DIR}/${DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}"
  done
}

for thermo in cmp fmt hmw ymp ypf; do
  declare -a problems
  case "${thermo}" in
    "cmp")
      problems=(acidmwb bo3bufs cahco3 co3aqui cristoba custbuf deionw fo2mineq h2so4p02 h2so4p1 hcl1 hclp02 hclp1 henleyph j13w j13wa j13wsf oxcalhem ph4hcl quenchfl rwpar rwssdiag rwtst sio2 sw2vxch swauto swbasw swcaarag swco2 swg1sx swgeo swmaj swmajd swmjdofe swpar swpharag swrdxcp swtst swv1sx swv1sxca swv2sx)
      ;;
    "fmt")
      problems=(c4pgwbN2 deadseaw f24vc7b3 f24vc7k4 f24vc7m gypnaclx swmajm)
      ;;
    "hmw")
      problems=(calnacl canshbr deadsea epsomite gypsum swlmahcl swmajp swphcl swv1sxk)
      ;;
    "ymp")
      problems=(acidmwb bo3bufs cahco3 co3aqui cristoba custbuf deionw fo2mineq h2so4p02 h2so4p1 hcl1 hclp02 hclp1 henleyph j13w j13wa j13wsf oxcalhem ph4hcl quenchfl rwpar rwssdiag rwtst sio2 sw2vxch swauto swbasw swcaarag swco2 swg1sx swgeo swmaj swmajd swmjdofe swpar swpharag swrdxcp swtst swv1sx swv1sxca swv2sx)
      ;;
    "ypf")
      problems=(Sylcacl2 arcmir arcsyl hal100 syl100 sylhal sylhal2)
      ;;
  esac

  for problem in "${problems[@]}"; do
    bats_test_function --description "Correct output (${thermo}/${problem})" -- check_output "${thermo}" "${problem}"
  done
done
