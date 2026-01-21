setup() {
  load 'test_helper/common-setup'
  _common_setup

  cp "${TEST_TMPDIR}/bin/eqpt" .
  cp "${TEST_TMPDIR}/bin/eq6" .
}

@test "Binary Exists" {
  assert_exists eq6
}

@test "Correctly replaces extension" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/crisqtz.6i" .
  ./eq6 cmp.d1 crisqtz.6i
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

@test "Actual extension is not relevant" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/crisqtz.6i" crisqtz.prob
  ./eq6 cmp.d1 crisqtz.prob
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

@test "Handles no extension" {
  cp "${BATS_TEST_DIRNAME}/data/eqpt/cmp.d1" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/cmp/crisqtz.6i" crisqtz
  ./eq6 cmp.d1 crisqtz
  assert_exists crisqtz.6o
  assert_exists crisqtz.6p
  assert_exists crisqtz.6tx
}

check_output() {
  local -r DIR="${1}"
  local -r PROBLEM="${2}"
  local -r SNAPSHOT_DIR="$(_snapshot_dir eq6)"
  shift 2

  cp "${BATS_TEST_DIRNAME}/data/eqpt/${DIR}.d0" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/${DIR}/${PROBLEM}.6i" .
  for ext in 6i 6o 6p; do
    cp "${BATS_TEST_DIRNAME}/data/eq6/${DIR}/${PROBLEM}.${ext}" "expected.${ext}"
  done

  run ./eqpt "${DIR}.d0"
  run ./eq6 "${DIR}.d1" "${PROBLEM}.6i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.6o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.6o

  # DGM: We are not going to check the 6o files for a couple of reasons:
  #        1. They can differ in the number of iterations taken. While that
  #           matters, it also kinda doesn't...
  #        2. We are checking the pickup, 6tx and csv, which have all of the
  #           the majority of the computed values of note.
  for ext in 6p 6tx; do
    assert_files_almost_same  "expected.${ext}" "${PROBLEM}.${ext}" "${@}"
  done

  for ext in 6o 6p 6tx; do
    assert_files_equal "${SNAPSHOT_DIR}/${DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}"
  done
}

check_snapshot_only() {
  local -r DIR="${1}"
  local -r PROBLEM="${2}"
  local -r SNAPSHOT_DIR="$(_snapshot_dir eq6)"
  shift 2

  cp "${BATS_TEST_DIRNAME}/data/eqpt/${DIR}.d0" .
  cp "${BATS_TEST_DIRNAME}/data/eq6/${DIR}/${PROBLEM}.6i" .
  for ext in 6i 6o 6p; do
    cp "${BATS_TEST_DIRNAME}/data/eq6/${DIR}/${PROBLEM}.${ext}" "expected.${ext}"
  done

  run ./eqpt "${DIR}.d0"
  run ./eq6 "${DIR}.d1" "${PROBLEM}.6i"

  perl -ni -e 'print unless /^\s*(Start|End|Run)\s+time/' ./*.6o
  perl -ni -e 'print unless /^\s*Run\s+[0-9]+/' ./*.6o

  for ext in 6o 6p 6tx; do
    assert_files_equal "${SNAPSHOT_DIR}/${DIR}/${PROBLEM}.${ext}" "${PROBLEM}.${ext}"
  done
}

@test "Correct output (cmp/crisqtz)" {
  case "$(arch)" in
    "arm64")
      check_output cmp crisqtz 1e-8
      ;;
    "x86_64"|"amd64")
      check_output cmp crisqtz 1e-14
      ;;
    *)
      check_output cmp crisqtz
      ;;
  esac
}
@test "Correct output (cmp/dedolo)" {
  case "$(arch)" in
    "arm64")
      check_output cmp dedolo 1e-7
      ;;
    *)
      check_output cmp dedolo
      ;;
  esac
}
@test "Correct output (cmp/heatqf)" {
  check_output cmp heatqf
}
@test "Correct output (cmp/heatsw)" {
  case "$(arch)" in
    "arm64")
      check_output cmp heatsw 1e-15
      ;;
    "x86_64"|"amd64")
      check_output cmp heatsw 1e-15
      ;;
    *)
      check_output cmp heatsw
      ;;
  esac
}
@test "Correct output (cmp/heatswfl)" {
  case "$(arch)" in
    "arm64")
      check_output cmp heatswfl 1e-5
      ;;
    "x86_64"|"amd64")
      check_output cmp heatswfl 1e-6
      ;;
    *)
      check_output cmp heatswfl
      ;;
  esac
}
@test "Correct output (cmp/j13wsf)" {
  check_output cmp j13wsf
}
@test "Correct output (cmp/j13wtitr)" {
  case "$(arch)" in
    "arm64")
      check_output cmp j13wtitr 1e-14
      ;;
    "x86_64"|"amd64")
      check_output cmp j13wtitr 1e-14
      ;;
    *)
      check_output cmp j13wtitr
      ;;
  esac
}
@test "Correct output (cmp/j13wtuff)" {
  case "$(arch)" in
    "arm64")
      check_output cmp j13wtuff 1e-9
      ;;
    "x86_64"|"amd64")
      check_output cmp j13wtuff 1e-9
      ;;
    *)
      check_output cmp j13wtuff
      ;;
  esac
}
@test "Correct output (cmp/methane)" {
  case "$(arch)" in
    "arm64")
      check_output cmp methane 1e-4
      ;;
    *)
      check_output cmp methane
      ;;
  esac
}
@test "Correct output (cmp/micro)" {
  case "$(arch)" in
    "arm64")
      check_output cmp micro 1e-7
      ;;
    *)
      check_output cmp micro
      ;;
  esac
}
@test "Correct output (cmp/microft)" {
  case "$(arch)" in
    "arm64")
      check_output cmp microft 1e-3
      ;;
    *)
      check_output cmp microft
      ;;
  esac
}
@test "Correct output (cmp/pptcal)" {
  case "$(arch)" in
    "arm64")
      check_output cmp pptcal 1e-15
      ;;
    *)
      check_output cmp pptcal
      ;;
  esac
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
  case "$(arch)" in
    "x86_64"|"amd64")
      check_output cmp pyrsw 1e-19
      ;;
    *)
      check_output cmp pyrsw
      ;;
  esac
}
@test "Correct output (cmp/rwssdiag)" {
  case "$(arch)" in
    "arm64")
      check_output cmp rwssdiag 1e-8
      ;;
    "x86_64"|"amd64")
      check_output cmp rwssdiag 1e-8
      ;;
    *)
      check_output cmp rwssdiag
      ;;
  esac
}
@test "Correct output (cmp/rwtitr)" {
  case "$(arch)" in
    "arm64")
      check_output cmp rwtitr 1e-15
      ;;
    "x86_64"|"amd64")
      check_output cmp rwtitr 1e-17
      ;;
    *)
      check_output cmp rwtitr
      ;;
  esac
}
@test "Correct output (cmp/swtitr)" {
  case "$(arch)" in
    "arm64")
      check_output cmp swtitr 1e-9
      ;;
    "x86_64"|"amd64")
      check_output cmp swtitr 1e-9
      ;;
    *)
      check_output cmp swtitr
      ;;
  esac
}
@test "Correct output (cmp/swxrca)" {
  case "$(arch)" in
    "arm64")
      skip "This problem is known to be broken on arm64"
      ;;
    "x86_64"|"amd64")
      check_output cmp swxrca 1e-6
      ;;
    *)
      check_output cmp swxrca
      ;;
  esac
}
@test "Check snapshot only (cmp/swxrca)" {
  check_snapshot_only cmp swxrca
}
@test "Correct output (cmp/swxrcaft)" {
  case "$(arch)" in
    "arm64")
      skip "This problem is known to be broken on arm64"
      ;;
    "x86_64"|"amd64")
      check_output cmp swxrcaft 1e-5
      ;;
    *)
      check_output cmp swxrcaft
      ;;
  esac
}
@test "Check snapshot only (cmp/swxrcaft)" {
  check_snapshot_only cmp swxrcaft
}
@test "Correct output (fmt/c4pgwbN2)" {
  case "$(arch)" in
    "arm64")
      check_output fmt c4pgwbN2 1e-11
      ;;
    "x86_64"|"amd64")
      check_output fmt c4pgwbN2 1e-11
      ;;
    *)
      check_output fmt c4pgwbN2
      ;;
  esac
}
@test "Correct output (fmt/f24vc7b3)" {
  case "$(arch)" in
    "arm64")
      check_output fmt f24vc7b3 1e-13
      ;;
    "x86_64"|"amd64")
      check_output fmt f24vc7b3 1e-13
      ;;
    *)
      check_output fmt f24vc7b3
      ;;
  esac
}
@test "Correct output (fmt/f24vc7k4)" {
  case "$(arch)" in
    "arm64")
      check_output fmt f24vc7k4 1e-11
      ;;
    "x86_64"|"amd64")
      check_output fmt f24vc7k4 1e-11
      ;;
    *)
      check_output fmt f24vc7k4
      ;;
  esac
}
@test "Correct output (fmt/f24vc7m)" {
  case "$(arch)" in
    "arm64")
      check_output fmt f24vc7m 1e-13
      ;;
    "x86_64"|"amd64")
      check_output fmt f24vc7m 1e-13
      ;;
    *)
      check_output fmt f24vc7m
      ;;
  esac
}
@test "Correct output (fmt/gypnaclx)" {
  case "$(arch)" in
    "arm64")
      check_output fmt gypnaclx 1e-13
      ;;
    *)
      check_output fmt gypnaclx
      ;;
  esac
}
@test "Correct output (hmw/calhal)" {
  case "$(arch)" in
    "arm64")
      check_output hmw calhal 1e-14
      ;;
    "x86_64"|"amd64")
      check_output hmw calhal 1e-14
      ;;
    *)
      check_output hmw calhal
      ;;
  esac
}
@test "Correct output (hmw/evapsw)" {
  case "$(arch)" in
    "arm64")
      check_output hmw evapsw 1e-6
      ;;
    "x86_64"|"amd64")
      check_output hmw evapsw 1e-7
      ;;
    *)
      check_output hmw evapsw
      ;;
  esac
}
@test "Correct output (hmw/evswgyha)" {
  case "$(arch)" in
    "arm64")
      check_output hmw evswgyha 1e-13
      ;;
    "x86_64"|"amd64")
      check_output hmw evswgyha 1e-12
      ;;
    *)
      check_output hmw evswgyha
      ;;
  esac
}
@test "Correct output (hmw/fwbrmix)" {
  case "$(arch)" in
    "arm64")
      check_output hmw fwbrmix 1e-13
      ;;
    "x86_64"|"amd64")
      check_output hmw fwbrmix 1e-14
      ;;
    *)
      check_output hmw fwbrmix
      ;;
  esac
}
@test "Correct output (hmw/gypanhy)" {
  case "$(arch)" in
    "arm64")
      check_output hmw gypanhy 1e-14
      ;;
    *)
      check_output hmw gypanhy
      ;;
  esac
}
@test "Correct output (hmw/mgso4)" {
  case "$(arch)" in
    "arm64")
      skip "This problem is known to be broken on arm64"
      ;;
    "x86_64"|"amd64")
      check_output hmw mgso4 1e-8
      ;;
    *)
      check_output hmw mgso4
      ;;
  esac
}
@test "Check snapshot only (hmw/mgso4)" {
  check_snapshot_only hmw mgso4
}
@test "Correct output (hmw/swv1sxk)" {
  case "$(arch)" in
    "arm64")
      check_output hmw swv1sxk 1e-15
      ;;
    *)
      check_output hmw swv1sxk
      ;;
  esac
}
@test "Correct output (ymp/crisqtz)" {
  case "$(arch)" in
    "arm64")
      check_output ymp crisqtz 1e-13
      ;;
    *)
      check_output ymp crisqtz
      ;;
  esac
}
@test "Correct output (ymp/dedolo)" {
  case "$(arch)" in
    "arm64")
      check_output ymp dedolo 1e-9
      ;;
    "x86_64"|"amd64")
      check_output ymp dedolo 1e-13
      ;;
    *)
      check_output ymp dedolo
      ;;
  esac
}
@test "Correct output (ymp/heatqf)" {
  case "$(arch)" in
    "arm64")
      check_output ymp heatqf 1e-18
      ;;
    *)
      check_output ymp heatqf
      ;;
  esac
}
@test "Correct output (ymp/heatsw)" {
  case "$(arch)" in
    "arm64")
      check_output ymp heatsw 1e-16
      ;;
    "x86_64"|"amd64")
      check_output ymp heatsw 1e-16
      ;;
    *)
      check_output ymp heatsw
      ;;
  esac
}
@test "Correct output (ymp/j13wsf)" {
  case "$(arch)" in
    "arm64")
      check_output ymp j13wsf 1e-6
      ;;
    "x86_64"|"amd64")
      check_output ymp j13wsf 1e-9
      ;;
    *)
      check_output ymp j13wsf
      ;;
  esac
}
@test "Correct output (ymp/j13wtitr)" {
  case "$(arch)" in
    "arm64")
      check_output ymp j13wtitr 1e-15
      ;;
    *)
      check_output ymp j13wtitr
      ;;
  esac
}
@test "Correct output (ymp/j13wtuff)" {
  case "$(arch)" in
    "arm64")
      check_output ymp j13wtuff 1e-6
      ;;
    "x86_64"|"amd64")
      check_output ymp j13wtuff 1e-6
      ;;
    *)
      check_output ymp j13wtuff
      ;;
  esac
}
@test "Correct output (ymp/methane)" {
  case "$(arch)" in
    "arm64")
      check_output ymp methane 1e-3
      ;;
    "x86_64"|"amd64")
      check_output ymp methane 1e-5
      ;;
    *)
      check_output ymp methane
      ;;
  esac
}
@test "Correct output (ymp/micro)" {
  case "$(arch)" in
    "arm64")
      check_output ymp micro 1e-13
      ;;
    "x86_64"|"amd64")
      check_output ymp micro 1e-13
      ;;
    *)
      check_output ymp micro
      ;;
  esac
}
@test "Correct output (ymp/microft)" {
  case "$(arch)" in
    "arm64")
      check_output ymp microft 1e-4
      ;;
    "x86_64"|"amd64")
      check_output ymp microft 1e-4
      ;;
    *)
      check_output ymp microft
      ;;
  esac
}
@test "Correct output (ymp/pptcal)" {
  case "$(arch)" in
    "arm64")
      check_output ymp pptcal 1e-15
      ;;
    *)
      check_output ymp pptcal
      ;;
  esac
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
  case "$(arch)" in
    "arm64")
      check_output ymp rwssdiag 1e-9
      ;;
    *)
      check_output ymp rwssdiag
      ;;
  esac
}
@test "Correct output (ymp/rwtitr)" {
  case "$(arch)" in
    "arm64")
      check_output ymp rwtitr 1e-14
      ;;
    *)
      check_output ymp rwtitr
      ;;
  esac
}
@test "Correct output (ymp/swtitr)" {
  case "$(arch)" in
    "arm64")
      check_output ymp swtitr 1e-13
      ;;
    "x86_64"|"amd64")
      check_output ymp swtitr 1e-13
      ;;
    *)
      check_output ymp swtitr
      ;;
  esac
}
@test "Correct output (ymp/swxrca)" {
  case "$(arch)" in
    "arm64")
      skip "This problem is known to be broken on arm64"
      ;;
    "x86_64"|"amd64")
      check_output ymp swxrca 1e-13
      ;;
    *)
      check_output ymp swxrca
      ;;
  esac
}
@test "Check snapshot only (ymp/swxrca)" {
  check_snapshot_only ymp swxrca
}
@test "Correct output (ymp/swxrcaft)" {
  case "$(arch)" in
    "arm64")
      skip "This problem is known to be broken on arm64"
      ;;
    "x86_64"|"amd64")
      check_output ymp swxrcaft 1e-8
      ;;
    *)
      check_output ymp swxrcaft
      ;;
  esac
}
@test "Check snapshot only (ymp/swxrcaft)" {
  check_snapshot_only ymp swxrcaft
}
@test "Correct output (ypf/calhal90)" {
  case "$(arch)" in
    "arm64")
      check_output ypf calhal90 1e-13
      ;;
    *)
      check_output ypf calhal90
      ;;
  esac
}
@test "Correct output (ypf/evapsw60)" {
  case "$(arch)" in
    "arm64")
      check_output ypf evapsw60 1e-11
      ;;
    "x86_64"|"amd64")
      check_output ypf evapsw60 1e-12
      ;;
    *)
      check_output ypf evapsw60
      ;;
  esac
}
