setup_file() {
  load '../../test_helper/common-setup'

  _create_snapshot_dir eq6
}

setup() {
  load '../../test_helper/common-setup'

  _snapshot_setup
}

snapshot() {
  cp "${BATS_TEST_DIRNAME}/../eqpt/${1}.d0" .
  cp "${BATS_TEST_DIRNAME}/../eq6/${1}/${2}.6i" .

  ./eqpt "${1}.d0"
  ./eq6 "${1}.d1" "${2}.6i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.6o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.6o

  mkdir -p "${SNAPSHOT_DIR}/${1}"
  for ext in 6o 6p 6tx; do
    cp "${2}.${ext}" "${SNAPSHOT_DIR}/${1}/${2}.${ext}"
  done
}

@test "Snapshot (cmp/crisqtz)" {
  snapshot cmp crisqtz
}
@test "Snapshot (cmp/dedolo)" {
  snapshot cmp dedolo
}
@test "Snapshot (cmp/heatqf)" {
  snapshot cmp heatqf
}
@test "Snapshot (cmp/heatsw)" {
  snapshot cmp heatsw
}
@test "Snapshot (cmp/heatswfl)" {
  snapshot cmp heatswfl
}
@test "Snapshot (cmp/j13wsf)" {
  snapshot cmp j13wsf
}
@test "Snapshot (cmp/j13wtitr)" {
  snapshot cmp j13wtitr
}
@test "Snapshot (cmp/j13wtuff)" {
  snapshot cmp j13wtuff
}
@test "Snapshot (cmp/methane)" {
  snapshot cmp methane
}
@test "Snapshot (cmp/micro)" {
  snapshot cmp micro
}
@test "Snapshot (cmp/microft)" {
  snapshot cmp microft
}
@test "Snapshot (cmp/pptcal)" {
  snapshot cmp pptcal
}
@test "Snapshot (cmp/pptmins)" {
  snapshot cmp pptmins
}
@test "Snapshot (cmp/pptqtz)" {
  snapshot cmp pptqtz
}
@test "Snapshot (cmp/pptqtza)" {
  snapshot cmp pptqtza
}
@test "Snapshot (cmp/pyrsw)" {
  snapshot cmp pyrsw
}
@test "Snapshot (cmp/rwssdiag)" {
  snapshot cmp rwssdiag
}
@test "Snapshot (cmp/rwtitr)" {
  snapshot cmp rwtitr
}
@test "Snapshot (cmp/swtitr)" {
  snapshot cmp swtitr
}
@test "Snapshot (cmp/swxrca)" {
  snapshot cmp swxrca
}
@test "Snapshot (cmp/swxrcaft)" {
  snapshot cmp swxrcaft
}
@test "Snapshot (fmt/c4pgwbN2)" {
  snapshot fmt c4pgwbN2
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
@test "Snapshot (hmw/calhal)" {
  snapshot hmw calhal
}
@test "Snapshot (hmw/evapsw)" {
  snapshot hmw evapsw
}
@test "Snapshot (hmw/evswgyha)" {
  snapshot hmw evswgyha
}
@test "Snapshot (hmw/fwbrmix)" {
  snapshot hmw fwbrmix
}
@test "Snapshot (hmw/gypanhy)" {
  snapshot hmw gypanhy
}
@test "Snapshot (hmw/mgso4)" {
  snapshot hmw mgso4
}
@test "Snapshot (hmw/swv1sxk)" {
  snapshot hmw swv1sxk
}
@test "Snapshot (ymp/crisqtz)" {
  snapshot ymp crisqtz
}
@test "Snapshot (ymp/dedolo)" {
  snapshot ymp dedolo
}
@test "Snapshot (ymp/heatqf)" {
  snapshot ymp heatqf
}
@test "Snapshot (ymp/heatsw)" {
  snapshot ymp heatsw
}
@test "Snapshot (ymp/j13wsf)" {
  snapshot ymp j13wsf
}
@test "Snapshot (ymp/j13wtitr)" {
  snapshot ymp j13wtitr
}
@test "Snapshot (ymp/j13wtuff)" {
  snapshot ymp j13wtuff
}
@test "Snapshot (ymp/methane)" {
  snapshot ymp methane
}
@test "Snapshot (ymp/micro)" {
  snapshot ymp micro
}
@test "Snapshot (ymp/microft)" {
  snapshot ymp microft
}
@test "Snapshot (ymp/pptcal)" {
  snapshot ymp pptcal
}
@test "Snapshot (ymp/pptmins)" {
  snapshot ymp pptmins
}
@test "Snapshot (ymp/pptqtz)" {
  snapshot ymp pptqtz
}
@test "Snapshot (ymp/pptqtza)" {
  snapshot ymp pptqtza
}
@test "Snapshot (ymp/pyrsw)" {
  snapshot ymp pyrsw
}
@test "Snapshot (ymp/rwssdiag)" {
  snapshot ymp rwssdiag
}
@test "Snapshot (ymp/rwtitr)" {
  snapshot ymp rwtitr
}
@test "Snapshot (ymp/swtitr)" {
  snapshot ymp swtitr
}
@test "Snapshot (ymp/swxrca)" {
  snapshot ymp swxrca
}
@test "Snapshot (ymp/swxrcaft)" {
  snapshot ymp swxrcaft
}
@test "Snapshot (ypf/calhal90)" {
  snapshot ypf calhal90
}
@test "Snapshot (ypf/evapsw60)" {
  snapshot ypf evapsw60
}
