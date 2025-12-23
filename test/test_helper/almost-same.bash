#!/usr/bin/env bash

function _compare() {
  local -r file1="${1}"
  local -r file2="${2}"

  local atol=${3}
  local rtol=${4}

  local atol=${atol:=${TEST_ATOL}}
  local rtol=${rtol:=${TEST_RTOL}}

  diff "${file1}" "${file2}" | "${PYTHON}" "${BATS_TEST_DIRNAME}/test_helper/almost-same.py" -a "${atol}" -r "${rtol}"
}

assert_files_almost_same() {
  local -r file1="${1}"
  local -r file2="${2}"

  run _compare "${@}"
  if [ "${status}" -ne 0 ]; then
    local -r rem="${BATSLIB_FILE_PATH_REM-}"
    local -r add="${BATSLIB_FILE_PATH_ADD-}"
    batslib_print_kv_single 6 'path' "${file1/${rem}/${add}}" 'path' "${file2/${rem}/${add}}" 'diff' "
${output}" \
      | batslib_decorate 'files are too different' \
      | fail
  fi
}
