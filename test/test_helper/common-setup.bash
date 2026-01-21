#!/usr/bin/env bash

function _common_setup() {
  load 'test_helper/bats-support/load'
  load 'test_helper/bats-assert/load'
  load 'test_helper/bats-file/load'
  load 'test_helper/almost-same'

  TEST_TMPDIR="$(cd "${BATS_TEST_DIRNAME}/.." >/dev/null 2>&1 && pwd)"
  export TEST_TMPDIR

  cd "${BATS_TEST_TMPDIR}" || return 1

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
