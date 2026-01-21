setup_file() {
  load '../../test_helper/common-setup'

  _create_snapshot_dir eq3nr
}

setup() {
  load '../../test_helper/common-setup'

  _snapshot_setup
}

snapshot() {
  cp "${BATS_TEST_DIRNAME}/../eqpt/${1}.d0" .
  cp "${BATS_TEST_DIRNAME}/../eq3nr/${1}/${2}.3i" .

  ./eqpt "${1}.d0"
  ./eq3nr "${1}.d1" "${2}.3i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.3o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.3o

  mkdir -p "${SNAPSHOT_DIR}/${1}"
  for ext in 3o 3p; do
    cp "${2}.${ext}" "${SNAPSHOT_DIR}/${1}/${2}.${ext}"
  done
}

@test "Snapshot (cmp/acidmwb)" {
  snapshot cmp acidmwb
}
@test "Snapshot (cmp/bo3bufs)" {
  snapshot cmp bo3bufs
}
@test "Snapshot (cmp/cahco3)" {
  snapshot cmp cahco3
}
@test "Snapshot (cmp/co3aqui)" {
  snapshot cmp co3aqui
}
@test "Snapshot (cmp/cristoba)" {
  snapshot cmp cristoba
}
@test "Snapshot (cmp/custbuf)" {
  snapshot cmp custbuf
}
@test "Snapshot (cmp/deionw)" {
  snapshot cmp deionw
}
@test "Snapshot (cmp/fo2mineq)" {
  snapshot cmp fo2mineq
}
@test "Snapshot (cmp/h2so4p02)" {
  snapshot cmp h2so4p02
}
@test "Snapshot (cmp/h2so4p1)" {
  snapshot cmp h2so4p1
}
@test "Snapshot (cmp/hcl1)" {
  snapshot cmp hcl1
}
@test "Snapshot (cmp/hclp02)" {
  snapshot cmp hclp02
}
@test "Snapshot (cmp/hclp1)" {
  snapshot cmp hclp1
}
@test "Snapshot (cmp/henleyph)" {
  snapshot cmp henleyph
}
@test "Snapshot (cmp/j13w)" {
  snapshot cmp j13w
}
@test "Snapshot (cmp/j13wa)" {
  snapshot cmp j13wa
}
@test "Snapshot (cmp/j13wsf)" {
  snapshot cmp j13wsf
}
@test "Snapshot (cmp/oxcalhem)" {
  snapshot cmp oxcalhem
}
@test "Snapshot (cmp/ph4hcl)" {
  snapshot cmp ph4hcl
}
@test "Snapshot (cmp/quenchfl)" {
  snapshot cmp quenchfl
}
@test "Snapshot (cmp/rwpar)" {
  snapshot cmp rwpar
}
@test "Snapshot (cmp/rwssdiag)" {
  snapshot cmp rwssdiag
}
@test "Snapshot (cmp/rwtst)" {
  snapshot cmp rwtst
}
@test "Snapshot (cmp/sio2)" {
  snapshot cmp sio2
}
@test "Snapshot (cmp/sw2vxch)" {
  snapshot cmp sw2vxch
}
@test "Snapshot (cmp/swauto)" {
  snapshot cmp swauto
}
@test "Snapshot (cmp/swbasw)" {
  snapshot cmp swbasw
}
@test "Snapshot (cmp/swcaarag)" {
  snapshot cmp swcaarag
}
@test "Snapshot (cmp/swco2)" {
  snapshot cmp swco2
}
@test "Snapshot (cmp/swg1sx)" {
  snapshot cmp swg1sx
}
@test "Snapshot (cmp/swgeo)" {
  snapshot cmp swgeo
}
@test "Snapshot (cmp/swmaj)" {
  snapshot cmp swmaj
}
@test "Snapshot (cmp/swmajd)" {
  snapshot cmp swmajd
}
@test "Snapshot (cmp/swmjdofe)" {
  snapshot cmp swmjdofe
}
@test "Snapshot (cmp/swpar)" {
  snapshot cmp swpar
}
@test "Snapshot (cmp/swpharag)" {
  snapshot cmp swpharag
}
@test "Snapshot (cmp/swrdxcp)" {
  snapshot cmp swrdxcp
}
@test "Snapshot (cmp/swtst)" {
  snapshot cmp swtst
}
@test "Snapshot (cmp/swv1sx)" {
  snapshot cmp swv1sx
}
@test "Snapshot (cmp/swv1sxca)" {
  snapshot cmp swv1sxca
}
@test "Snapshot (cmp/swv2sx)" {
  snapshot cmp swv2sx
}
@test "Snapshot (fmt/c4pgwbN2)" {
  snapshot fmt c4pgwbN2
}
@test "Snapshot (fmt/deadseaw)" {
  snapshot fmt deadseaw
}
@test "Snapshot (fmt/f24vc7b3)" {
  snapshot fmt f24vc7b3
}
@test "Snapshot (fmt/f24vc7k4)" {
  snapshot fmt f24vc7k4
}
@test "Snapshot (fmt/f24vc7m)" {
  snapshot fmt f24vc7m
}
@test "Snapshot (fmt/gypnaclx)" {
  snapshot fmt gypnaclx
}
@test "Snapshot (fmt/swmajm)" {
  snapshot fmt swmajm
}
@test "Snapshot (hmw/calnacl)" {
  snapshot hmw calnacl
}
@test "Snapshot (hmw/canshbr)" {
  snapshot hmw canshbr
}
@test "Snapshot (hmw/deadsea)" {
  snapshot hmw deadsea
}
@test "Snapshot (hmw/epsomite)" {
  snapshot hmw epsomite
}
@test "Snapshot (hmw/gypsum)" {
  snapshot hmw gypsum
}
@test "Snapshot (hmw/swlmahcl)" {
  snapshot hmw swlmahcl
}
@test "Snapshot (hmw/swmajp)" {
  snapshot hmw swmajp
}
@test "Snapshot (hmw/swphcl)" {
  snapshot hmw swphcl
}
@test "Snapshot (hmw/swv1sxk)" {
  snapshot hmw swv1sxk
}
@test "Snapshot (ymp/acidmwb)" {
  snapshot ymp acidmwb
}
@test "Snapshot (ymp/bo3bufs)" {
  snapshot ymp bo3bufs
}
@test "Snapshot (ymp/cahco3)" {
  snapshot ymp cahco3
}
@test "Snapshot (ymp/co3aqui)" {
  snapshot ymp co3aqui
}
@test "Snapshot (ymp/cristoba)" {
  snapshot ymp cristoba
}
@test "Snapshot (ymp/custbuf)" {
  snapshot ymp custbuf
}
@test "Snapshot (ymp/deionw)" {
  snapshot ymp deionw
}
@test "Snapshot (ymp/fo2mineq)" {
  snapshot ymp fo2mineq
}
@test "Snapshot (ymp/h2so4p02)" {
  snapshot ymp h2so4p02
}
@test "Snapshot (ymp/h2so4p1)" {
  snapshot ymp h2so4p1
}
@test "Snapshot (ymp/hcl1)" {
  snapshot ymp hcl1
}
@test "Snapshot (ymp/hclp02)" {
  snapshot ymp hclp02
}
@test "Snapshot (ymp/hclp1)" {
  snapshot ymp hclp1
}
@test "Snapshot (ymp/henleyph)" {
  snapshot ymp henleyph
}
@test "Snapshot (ymp/j13w)" {
  snapshot ymp j13w
}
@test "Snapshot (ymp/j13wa)" {
  snapshot ymp j13wa
}
@test "Snapshot (ymp/j13wsf)" {
  snapshot ymp j13wsf
}
@test "Snapshot (ymp/oxcalhem)" {
  snapshot ymp oxcalhem
}
@test "Snapshot (ymp/ph4hcl)" {
  snapshot ymp ph4hcl
}
@test "Snapshot (ymp/quenchfl)" {
  snapshot ymp quenchfl
}
@test "Snapshot (ymp/rwpar)" {
  snapshot ymp rwpar
}
@test "Snapshot (ymp/rwssdiag)" {
  snapshot ymp rwssdiag
}
@test "Snapshot (ymp/rwtst)" {
  snapshot ymp rwtst
}
@test "Snapshot (ymp/sio2)" {
  snapshot ymp sio2
}
@test "Snapshot (ymp/sw2vxch)" {
  snapshot ymp sw2vxch
}
@test "Snapshot (ymp/swauto)" {
  snapshot ymp swauto
}
@test "Snapshot (ymp/swbasw)" {
  snapshot ymp swbasw
}
@test "Snapshot (ymp/swcaarag)" {
  snapshot ymp swcaarag
}
@test "Snapshot (ymp/swco2)" {
  snapshot ymp swco2
}
@test "Snapshot (ymp/swg1sx)" {
  snapshot ymp swg1sx
}
@test "Snapshot (ymp/swgeo)" {
  snapshot ymp swgeo
}
@test "Snapshot (ymp/swmaj)" {
  snapshot ymp swmaj
}
@test "Snapshot (ymp/swmajd)" {
  snapshot ymp swmajd
}
@test "Snapshot (ymp/swmjdofe)" {
  snapshot ymp swmjdofe
}
@test "Snapshot (ymp/swpar)" {
  snapshot ymp swpar
}
@test "Snapshot (ymp/swpharag)" {
  snapshot ymp swpharag
}
@test "Snapshot (ymp/swrdxcp)" {
  snapshot ymp swrdxcp
}
@test "Snapshot (ymp/swtst)" {
  snapshot ymp swtst
}
@test "Snapshot (ymp/swv1sx)" {
  snapshot ymp swv1sx
}
@test "Snapshot (ymp/swv1sxca)" {
  snapshot ymp swv1sxca
}
@test "Snapshot (ymp/swv2sx)" {
  snapshot ymp swv2sx
}
@test "Snapshot (ypf/Sylcacl2)" {
  snapshot ypf Sylcacl2
}
@test "Snapshot (ypf/arcmir)" {
  snapshot ypf arcmir
}
@test "Snapshot (ypf/arcsyl)" {
  snapshot ypf arcsyl
}
@test "Snapshot (ypf/hal100)" {
  snapshot ypf hal100
}
@test "Snapshot (ypf/syl100)" {
  snapshot ypf syl100
}
@test "Snapshot (ypf/sylhal)" {
  snapshot ypf sylhal
}
@test "Snapshot (ypf/sylhal2)" {
  snapshot ypf sylhal2
}
