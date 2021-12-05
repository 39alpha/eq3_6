#!/usr/bin/env bash

TIMING_REGEX='((Start|End|Run)[[:space:]]+time|Run[[:space:]]+[[:digit:]]{2})'

failures=0
tests=0
failedtests=()

almost_same="$PYTHON $TEST_ROOT/almost_same.py"

function compare() {
    if [ $# -ge 3 ]; then
        awk "!/$3/{print}" $1 > out.a
        awk "!/$3/{print}" $2 > out.b
        diff out.a out.b | $almost_same
    else
        diff $1 $2 | $almost_same
    fi
}

function assert_same_contents() {
    if [ $# -ge 3 ]; then
        OUTPUT=$(compare $1 $2 $TIMING_REGEX)
    else
        OUTPUT=$(compare $@)
    fi

    if [ $? -ne 0 ]; then
        echo "FAILED: checking \"$1\" against \"$2\"" >&2
        echo -e "$OUTPUT\n" >&2
        failures=$(($failures + 1))
        failedtests+=($groundtruth)
    fi
}

function summarize() {
    echo "SUMMARY"
    if (($failures != 0)); then
        for tst in "${failedtests[@]}"; do
            if [ -f "$tst" ]; then
                echo "  FAILED: $(realpath --relative-to=.. $tst)"
            else
                echo "  FAILED: $tst"
            fi
        done
        echo ''
    fi
    echo "  $failures of $tests tests failed"

    exit $failures
}
