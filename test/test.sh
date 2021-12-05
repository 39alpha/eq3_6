#!/usr/bin/env bash

UNITS=(eqpt eq3nr eq6)

failed=0
for unit in "${UNITS[@]}"; do
    pushd $unit >/dev/null 2>/dev/null
    if ! bash test.sh; then
        failed=$((failed + 1))
    fi
    popd >/dev/null 2>/dev/null
done

exit $failed
