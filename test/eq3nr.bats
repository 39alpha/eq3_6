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

  run which eq3nr
  assert_output "${BATS_TEST_TMPDIR}/bin/eq3nr"

  run pwd
  assert_output "${BATS_TEST_TMPDIR}"
}

@test "Binary Exists" {
  is_setup
}

@test "Correctly replaces extension" {
  is_setup

  cp "${EQ36_TEST_DATA_DIR}/eqpt/cmp.d1" .
  cp "${EQ36_TEST_DATA_DIR}/eq3nr/cmp/acidmwb.3i" .
  eq3nr cmp.d1 acidmwb.3i
  assert_exists acidmwb.3o
  assert_exists acidmwb.3p
}

@test "Actual extension is not relevant" {
  is_setup

  cp "${EQ36_TEST_DATA_DIR}/eqpt/cmp.d1" .
  cp "${EQ36_TEST_DATA_DIR}/eq3nr/cmp/acidmwb.3i" acidmwb.prob
  eq3nr cmp.d1 acidmwb.prob
  assert_exists acidmwb.3o
  assert_exists acidmwb.3p
}

@test "Handles no extension" {
  is_setup

  cp "${EQ36_TEST_DATA_DIR}/eqpt/cmp.d1" .
  cp "${EQ36_TEST_DATA_DIR}/eq3nr/cmp/acidmwb.3i" acidmwb
  eq3nr cmp.d1 acidmwb
  assert_exists acidmwb.3o
  assert_exists acidmwb.3p
}

check_output() {
  is_setup

  local -r DIR="${1}"
  local -r PROBLEM="${2}"
  shift 2

  cp "${EQ36_TEST_DATA_DIR}/eqpt/${DIR}.d0" .
  cp "${EQ36_TEST_DATA_DIR}/eq3nr/${DIR}/${PROBLEM}.3i" .
  for ext in 3o 3p; do
    cp "${EQ36_TEST_DATA_DIR}/eq3nr/${DIR}/${PROBLEM}.${ext}" "expected.${ext}"
  done

  run eqpt "${DIR}.d0"
  run eq3nr "${DIR}.d1" "${PROBLEM}.3i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.3o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.3o

  assert_files_almost_same "expected.3p" "${PROBLEM}.3p" 1e-12 1e-6
  assert_files_almost_same "expected.3o" "${PROBLEM}.3o" 1e-12 1e-6
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
