#!/usr/bin/env bash

export PROJECT_ROOT=$(realpath $(pwd)/..)
export TEST_ROOT=$(pwd)

if ! which python3 >/dev/null; then
    if ! which python >/dev/null; then
        echo "ERROR: looks like we can't find your python executable. Is it on your path" >&2
        exit -1
    else
        export PYTHON=$(which python)
    fi
else
    export PYTHON=$(which python3)
fi

UNITS=(eq3nr)

failed=0
for unit in "${UNITS[@]}"; do
    pushd $unit >/dev/null 2>&1
    if ! bash test.sh; then
        failed=$((failed + 1))
    fi
    popd >/dev/null 2>&1
done

exit $failed
