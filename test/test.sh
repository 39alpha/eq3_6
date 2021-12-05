#!/usr/bin/env bash

export PROJECT_ROOT=$(realpath $(pwd)/..)
export TEST_ROOT=$(pwd)

UNITS=(eqpt eq3nr eq6)

failed=0
for unit in "${UNITS[@]}"; do
    pushd $unit >/dev/null 2>&1
    if ! bash test.sh; then
        failed=$((failed + 1))
    fi
    popd >/dev/null 2>&1
done

exit $failed
