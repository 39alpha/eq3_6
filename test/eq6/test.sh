#!/usr/bin/env bash

EXPECTED_FILES=(output pickup)
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

EQ6=$(realpath -m ../../bin/eq6)
if ! [ -f $EQ6 ]; then
    echo "eq6 binary not found at $EQ6; perhaps you should run make?"
    exit -1
fi

echo "Running tests for $EQ6"
rm -rf tmp
for dir in data/*; do
    dataset=$(basename $dir)
    for sixi in $dir/*.6i; do
        problem=$(basename $sixi .6i)
        mkdir tmp
        cp $dir/data1.$dataset tmp
        cp $sixi tmp
        cd tmp
        $EQ6 data1.$dataset $problem.6i >/dev/null 2>&1
        for file in "${EXPECTED_FILES[@]}"; do
            if [ $file == "output" ]; then
                groundtruth=$(realpath ../$dir/$problem.6o)
                OUTPUT=$(compare $groundtruth $file $TIMING_REGEX)
                if [ $? -ne 0 ]; then
                    echo "FAILED: checking \"$file\" against \"$groundtruth\" ignoring dates and times"
                    echo -e "$OUTPUT\n"
                    failures=$(($failures + 1))
                    failedtests+=($groundtruth)
                fi
            elif [ $file == "pickup" ]; then
                groundtruth=$(realpath ../$dir/$problem.6p)
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
