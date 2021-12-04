#!/usr/bin/env bash

UNITS=(eqpt)

failed=0
for unit in "${UNITS[@]}"; do
    pushd eqpt >/dev/null 2>/dev/null
    if ! bash test.sh; then
        failed=$((failed + 1))
    fi
    popd >/dev/null 2>/dev/null
done

exit $failed
