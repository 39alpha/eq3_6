#!/usr/bin/env bash

function find_project_root() {
  local root="${BATS_TEST_DIRNAME}"
  while [ ! -f "${root}/meson.build" ]; do
    root="${root}/.."
  done
  realpath "${root}"
}

function _suite_setup() {
  EQ36_ROOT_DIR=$(find_project_root)
  export EQ36_ROOT_DIR

  EQ36_TEST_DIR="${EQ36_ROOT_DIR}/test"
  export EQ36_TEST_DIR

  EQ36_TEST_DATA_DIR="${EQ36_TEST_DIR}/data"
  export EQ36_TEST_DATA_DIR

  EQ36_TEST_BUILD_DIR="${EQ36_ROOT_DIR}/build/test"
  export EQ36_TEST_BUILD_DIR

  if [ -n "${EQ36_SNAPSHOTS:-}" ]; then
    EQ36_TEST_SNAPSHOTS_DIR="${EQ36_TEST_DATA_DIR}/snapshots/${EQ36_SNAPSHOTS}"
    export EQ36_TEST_SNAPSHOTS_DIR
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
  load "${EQ36_TEST_DIR}/test_helper/bats-support/load"
  load "${EQ36_TEST_DIR}/test_helper/bats-assert/load"
  load "${EQ36_TEST_DIR}/test_helper/bats-file/load"
  load "${EQ36_TEST_DIR}/test_helper/almost-same"

  mkdir -p "${BATS_TEST_TMPDIR}/bin"
  cp "${EQ36_TEST_BUILD_DIR}/eqpt" "${BATS_TEST_TMPDIR}/bin/eqpt"
  cp "${EQ36_TEST_BUILD_DIR}/eq3nr" "${BATS_TEST_TMPDIR}/bin/eq3nr"
  cp "${EQ36_TEST_BUILD_DIR}/eq6" "${BATS_TEST_TMPDIR}/bin/eq6"

  PATH="${BATS_TEST_TMPDIR}/bin:${PATH}"
  export PATH

  cd "${BATS_TEST_TMPDIR}" || return 2
}
