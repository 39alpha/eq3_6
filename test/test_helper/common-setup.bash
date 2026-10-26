#!/usr/bin/env bash

function _cli_build_dir() {
  if [ -n "${EQ36_TEST_BUILD_DIR:-}" ]; then
    echo "${EQ36_TEST_BUILD_DIR}"
    return
  fi
  local root
  root="$(cd "${BATS_TEST_DIRNAME}/../.." >/dev/null 2>&1 && pwd)"
  echo "${root}/build/test"
}

function _suite_setup() {
  local build root
  build="$(_cli_build_dir)"

  if [ -z "${EQ36_TEST_BUILD_DIR:-}" ]; then
    root="$(cd "${BATS_TEST_DIRNAME}/../.." >/dev/null 2>&1 && pwd)"
    mkdir -p "${root}/build"
    exec 9>"${root}/build/.cli-test.lock"
    flock 9
    [ -d "${build}" ] || meson setup "${build}" --buildtype=debug >/dev/null
    meson compile -C "${build}" >/dev/null
    flock -u 9
    exec 9>&-
  fi

  if ! which python3 >/dev/null; then
    if ! which python >/dev/null; then
      exit 1
    fi
    PYTHON=$(which python)
  else
    PYTHON=$(which python3)
  fi
  export PYTHON

  TEST_RTOL=${TEST_RTOL:=0}
  export TEST_RTOL

  TEST_ATOL=${TEST_ATOL:=0}
  export TEST_ATOL

  MALLOC_PERTURB_=${MALLOC_PERTURB_:=42}
  export MALLOC_PERTURB_
}

function _common_setup() {
  load 'test_helper/bats-support/load'
  load 'test_helper/bats-assert/load'
  load 'test_helper/bats-file/load'
  load 'test_helper/almost-same'

  TEST_TMPDIR="$(cd "${BATS_TEST_DIRNAME}/.." >/dev/null 2>&1 && pwd)"
  export TEST_TMPDIR

  local build
  build="$(_cli_build_dir)"
  mkdir -p "${BATS_TEST_TMPDIR}/bin"
  cp "${build}/eqpt" "${BATS_TEST_TMPDIR}/bin/eqpt"
  cp "${build}/eq3nr" "${BATS_TEST_TMPDIR}/bin/eq3nr"
  cp "${build}/eq6" "${BATS_TEST_TMPDIR}/bin/eq6"

  PATH="${BATS_TEST_TMPDIR}/bin:${PATH}"
  export PATH

  cd "${BATS_TEST_TMPDIR}" || return 2
}

function _snapshot_dir() {
  case "$(arch)" in
    "x86_64")
      local -r ARCH="amd64"
      ;;
    *)
      local -r ARCH="$(arch)"
      ;;
  esac

  if [[ -d "${BATS_TEST_DIRNAME}/data" ]]; then
    echo "${BATS_TEST_DIRNAME}/data/snapshots/${ARCH}/${*}"
  else
    echo "${BATS_TEST_DIRNAME}/${ARCH}/${*}"
  fi
}

function _create_snapshot_dir() {
  SNAPSHOT_DIR="$(_snapshot_dir "${@}")"

  rm -rf "${SNAPSHOT_DIR}"
  mkdir -p "${SNAPSHOT_DIR}"

  export SNAPSHOT_DIR
}

function _snapshot_setup() {
  load '../../test_helper/bats-support/load'
  load '../../test_helper/bats-assert/load'
  load '../../test_helper/bats-file/load'

  TEST_TMPDIR="$(cd "${BATS_TEST_DIRNAME}/../../.." >/dev/null 2>&1 && pwd)"
  export TEST_TMPDIR

  cd "${BATS_TEST_TMPDIR}" || return 1

  cp "${TEST_TMPDIR}/bin/eqpt" .
  cp "${TEST_TMPDIR}/bin/eq3nr" .
  cp "${TEST_TMPDIR}/bin/eq6" .
}
