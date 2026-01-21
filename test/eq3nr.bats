setup() {
  load 'test_helper/common-setup'
  _common_setup

  cp "${TEST_TMPDIR}/bin/eqpt" .
  cp "${TEST_TMPDIR}/bin/eq3nr" .
}

@test "Binary Exists" {
  assert_exists eq3nr
}

@test "Correctly replaces extension" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq3nr/cmp/acidmwb.3i" .
  ./eq3nr cmp.d1 acidmwb.3i
  assert_exists acidmwb.3o
  assert_exists acidmwb.3p
}

@test "Actual extension is not relevant" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq3nr/cmp/acidmwb.3i" acidmwb.prob
  ./eq3nr cmp.d1 acidmwb.prob
  assert_exists acidmwb.3o
  assert_exists acidmwb.3p
}

@test "Handles no extension" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq3nr/cmp/acidmwb.3i" acidmwb
  ./eq3nr cmp.d1 acidmwb
  assert_exists acidmwb.3o
  assert_exists acidmwb.3p
}

check_output() {
  local -r DIR="${1}"
  local -r PROBLEM="${2}"
  local -r SNAPSHOT_DIR="$(_snapshot_dir eq3nr)"
  shift 2

  cp "${BATS_TEST_DIRNAME}/data/eqpt/${DIR}.d0" .
  cp "${BATS_TEST_DIRNAME}/data/eq3nr/${DIR}/${PROBLEM}.3i" .
  for ext in 3o 3p; do
    cp "${BATS_TEST_DIRNAME}/data/eq3nr/${DIR}/${PROBLEM}.${ext}" "expected.${ext}"
  done

  run ./eqpt "${DIR}.d0"
  run ./eq3nr "${DIR}.d1" "${PROBLEM}.3i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.3o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.3o

  for ext in 3o 3p; do
    assert_files_almost_same "expected.${ext}" "${PROBLEM}.${ext}" "${@}"
  done

  for ext in 3o 3p; do
    assert_files_equal "${SNAPSHOT_DIR}/${DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}"
  done
}

@test "Correct output (cmp/acidmwb)" {
  check_output cmp acidmwb
}
@test "Correct output (cmp/bo3bufs)" {
  check_output cmp bo3bufs
}
@test "Correct output (cmp/cahco3)" {
  check_output cmp cahco3
}
@test "Correct output (cmp/co3aqui)" {
  check_output cmp co3aqui
}
@test "Correct output (cmp/cristoba)" {
  check_output cmp cristoba
}
@test "Correct output (cmp/custbuf)" {
  check_output cmp custbuf
}
@test "Correct output (cmp/deionw)" {
  check_output cmp deionw
}
@test "Correct output (cmp/fo2mineq)" {
  check_output cmp fo2mineq
}
@test "Correct output (cmp/h2so4p02)" {
  check_output cmp h2so4p02
}
@test "Correct output (cmp/h2so4p1)" {
  check_output cmp h2so4p1
}
@test "Correct output (cmp/hcl1)" {
  check_output cmp hcl1
}
@test "Correct output (cmp/hclp02)" {
  check_output cmp hclp02
}
@test "Correct output (cmp/hclp1)" {
  check_output cmp hclp1
}
@test "Correct output (cmp/henleyph)" {
  check_output cmp henleyph
}
@test "Correct output (cmp/j13w)" {
  check_output cmp j13w
}
@test "Correct output (cmp/j13wa)" {
  check_output cmp j13wa
}
@test "Correct output (cmp/j13wsf)" {
  check_output cmp j13wsf
}
@test "Correct output (cmp/oxcalhem)" {
  check_output cmp oxcalhem
}
@test "Correct output (cmp/ph4hcl)" {
  check_output cmp ph4hcl
}
@test "Correct output (cmp/quenchfl)" {
  check_output cmp quenchfl
}
@test "Correct output (cmp/rwpar)" {
  check_output cmp rwpar
}
@test "Correct output (cmp/rwssdiag)" {
  check_output cmp rwssdiag
}
@test "Correct output (cmp/rwtst)" {
  check_output cmp rwtst
}
@test "Correct output (cmp/sio2)" {
  check_output cmp sio2
}
@test "Correct output (cmp/sw2vxch)" {
  check_output cmp sw2vxch
}
@test "Correct output (cmp/swauto)" {
  check_output cmp swauto
}
@test "Correct output (cmp/swbasw)" {
  check_output cmp swbasw
}
@test "Correct output (cmp/swcaarag)" {
  check_output cmp swcaarag
}
@test "Correct output (cmp/swco2)" {
  check_output cmp swco2
}
@test "Correct output (cmp/swg1sx)" {
  check_output cmp swg1sx
}
@test "Correct output (cmp/swgeo)" {
  check_output cmp swgeo
}
@test "Correct output (cmp/swmaj)" {
  check_output cmp swmaj
}
@test "Correct output (cmp/swmajd)" {
  check_output cmp swmajd
}
@test "Correct output (cmp/swmjdofe)" {
  check_output cmp swmjdofe
}
@test "Correct output (cmp/swpar)" {
  check_output cmp swpar
}
@test "Correct output (cmp/swpharag)" {
  check_output cmp swpharag
}
@test "Correct output (cmp/swrdxcp)" {
  check_output cmp swrdxcp
}
@test "Correct output (cmp/swtst)" {
  check_output cmp swtst
}
@test "Correct output (cmp/swv1sx)" {
  check_output cmp swv1sx
}
@test "Correct output (cmp/swv1sxca)" {
  check_output cmp swv1sxca
}
@test "Correct output (cmp/swv2sx)" {
  check_output cmp swv2sx
}
@test "Correct output (fmt/c4pgwbN2)" {
  check_output fmt c4pgwbN2 1e-8
}
@test "Correct output (fmt/deadseaw)" {
  check_output fmt deadseaw
}
@test "Correct output (fmt/f24vc7b3)" {
  check_output fmt f24vc7b3
}
@test "Correct output (fmt/f24vc7k4)" {
  check_output fmt f24vc7k4
}
@test "Correct output (fmt/f24vc7m)" {
  check_output fmt f24vc7m
}
@test "Correct output (fmt/gypnaclx)" {
  check_output fmt gypnaclx
}
@test "Correct output (fmt/swmajm)" {
  check_output fmt swmajm
}
@test "Correct output (hmw/calnacl)" {
  check_output hmw calnacl
}
@test "Correct output (hmw/canshbr)" {
  check_output hmw canshbr
}
@test "Correct output (hmw/deadsea)" {
  check_output hmw deadsea
}
@test "Correct output (hmw/epsomite)" {
  check_output hmw epsomite
}
@test "Correct output (hmw/gypsum)" {
  check_output hmw gypsum
}
@test "Correct output (hmw/swlmahcl)" {
  check_output hmw swlmahcl
}
@test "Correct output (hmw/swmajp)" {
  check_output hmw swmajp
}
@test "Correct output (hmw/swphcl)" {
  check_output hmw swphcl
}
@test "Correct output (hmw/swv1sxk)" {
  check_output hmw swv1sxk
}
@test "Correct output (ymp/acidmwb)" {
  check_output ymp acidmwb
}
@test "Correct output (ymp/bo3bufs)" {
  check_output ymp bo3bufs
}
@test "Correct output (ymp/cahco3)" {
  check_output ymp cahco3
}
@test "Correct output (ymp/co3aqui)" {
  check_output ymp co3aqui
}
@test "Correct output (ymp/cristoba)" {
  check_output ymp cristoba
}
@test "Correct output (ymp/custbuf)" {
  check_output ymp custbuf
}
@test "Correct output (ymp/deionw)" {
  check_output ymp deionw
}
@test "Correct output (ymp/fo2mineq)" {
  check_output ymp fo2mineq
}
@test "Correct output (ymp/h2so4p02)" {
  check_output ymp h2so4p02
}
@test "Correct output (ymp/h2so4p1)" {
  check_output ymp h2so4p1
}
@test "Correct output (ymp/hcl1)" {
  check_output ymp hcl1
}
@test "Correct output (ymp/hclp02)" {
  check_output ymp hclp02
}
@test "Correct output (ymp/hclp1)" {
  check_output ymp hclp1
}
@test "Correct output (ymp/henleyph)" {
  check_output ymp henleyph
}
@test "Correct output (ymp/j13w)" {
  check_output ymp j13w
}
@test "Correct output (ymp/j13wa)" {
  check_output ymp j13wa
}
@test "Correct output (ymp/j13wsf)" {
  check_output ymp j13wsf
}
@test "Correct output (ymp/oxcalhem)" {
  check_output ymp oxcalhem
}
@test "Correct output (ymp/ph4hcl)" {
  check_output ymp ph4hcl
}
@test "Correct output (ymp/quenchfl)" {
  check_output ymp quenchfl
}
@test "Correct output (ymp/rwpar)" {
  check_output ymp rwpar
}
@test "Correct output (ymp/rwssdiag)" {
  check_output ymp rwssdiag
}
@test "Correct output (ymp/rwtst)" {
  check_output ymp rwtst
}
@test "Correct output (ymp/sio2)" {
  check_output ymp sio2
}
@test "Correct output (ymp/sw2vxch)" {
  check_output ymp sw2vxch
}
@test "Correct output (ymp/swauto)" {
  check_output ymp swauto
}
@test "Correct output (ymp/swbasw)" {
  check_output ymp swbasw
}
@test "Correct output (ymp/swcaarag)" {
  check_output ymp swcaarag
}
@test "Correct output (ymp/swco2)" {
  check_output ymp swco2
}
@test "Correct output (ymp/swg1sx)" {
  check_output ymp swg1sx
}
@test "Correct output (ymp/swgeo)" {
  check_output ymp swgeo
}
@test "Correct output (ymp/swmaj)" {
  check_output ymp swmaj
}
@test "Correct output (ymp/swmajd)" {
  check_output ymp swmajd
}
@test "Correct output (ymp/swmjdofe)" {
  check_output ymp swmjdofe
}
@test "Correct output (ymp/swpar)" {
  check_output ymp swpar
}
@test "Correct output (ymp/swpharag)" {
  check_output ymp swpharag
}
@test "Correct output (ymp/swrdxcp)" {
  check_output ymp swrdxcp
}
@test "Correct output (ymp/swtst)" {
  check_output ymp swtst
}
@test "Correct output (ymp/swv1sx)" {
  check_output ymp swv1sx
}
@test "Correct output (ymp/swv1sxca)" {
  check_output ymp swv1sxca
}
@test "Correct output (ymp/swv2sx)" {
  check_output ymp swv2sx
}
@test "Correct output (ypf/Sylcacl2)" {
  check_output ypf Sylcacl2
}
@test "Correct output (ypf/arcmir)" {
  check_output ypf arcmir 1e-15
}
@test "Correct output (ypf/arcsyl)" {
  check_output ypf arcsyl
}
@test "Correct output (ypf/hal100)" {
  check_output ypf hal100
}
@test "Correct output (ypf/syl100)" {
  check_output ypf syl100
}
@test "Correct output (ypf/sylhal)" {
  check_output ypf sylhal
}
@test "Correct output (ypf/sylhal2)" {
  check_output ypf sylhal2
}
