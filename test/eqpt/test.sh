#!/usr/bin/env bash

EXPECTED_FILES=(output data1 data1f slist)
TIMING_REGEX='((Start|End|Run)[[:space:]]+time|Run[[:space:]]+[[:digit:]]{2})'

failures=0
tests=0
failedtests=()

function compare() {
    if [ $# -eq 3 ]; then
        awk "!/$3/{print}" $1 > out.a
        awk "!/$3/{print}" $2 > out.b
        diff out.a out.b
    else
        diff $1 $2
    fi
}

EQPT=$(realpath -m ../../bin/eqpt)
if ! [ -f $EQPT ]; then
    echo "eqpt binary not found at $EQPT; perhaps you should run make?"
    exit -1
fi

echo "Running tests for $EQPT"
rm -rf tmp
for dir in data/*; do
    dataset=$(basename $dir)
    mkdir tmp
    cp $dir/data0.$dataset tmp
    cd tmp
    $EQPT data0.$dataset >/dev/null 2>&1
    for file in "${EXPECTED_FILES[@]}"; do
        groundtruth=$(realpath ../$dir/$file.$dataset)
        if [ $file == "output" ]; then
            OUTPUT=$(compare $groundtruth $file $TIMING_REGEX)
            if [ $? -ne 0 ]; then
                echo "FAILED: checking \"$file\" against \"$groundtruth\" ignoring dates and times"
                echo -e "$OUTPUT\n"
                failures=$(($failures + 1))
                failedtests+=($groundtruth)
            fi
        else
            OUTPUT=$(compare $groundtruth $file)
            if [ $? -ne 0 ]; then
                echo "FAILED: checking \"$file\" against \"$groundtruth\""
                echo -e "$OUTPUT\n"
                failures=$((failures + 1))
                failedtests+=($groundtruth)
            fi
        fi
        tests=$((tests + 1))
    done
    cd ..
    rm -rf tmp
done

echo "SUMMARY"
if (($failures != 0)); then
    for file in "${failedtests[@]}"; do
        echo "  FAILED $(realpath --relative-to=.. $file)"
    done
    echo ''
fi
echo "  $failures of $tests tests failed"

exit $failures
