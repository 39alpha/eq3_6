#!/usr/bin/env bash

source $TEST_ROOT/util.sh

EXTENSIONS=(3o 3p)

EQ3NR=$(realpath -m $PROJECT_ROOT/bin/eq3nr)
if ! [ -f $EQ3NR ]; then
    echo "eq3nr binary not found at $EQ3NR; perhaps you should run make?"
    exit -1
fi

echo "  [TEST] Correctly names files"
for name in acidmwb.3i acidwmb; do
    rm -rf tmp
    mkdir tmp
    cd tmp

    cp ../data/cmp/cmp.d1 .
    cp ../data/cmp/acidmwb.3i $name

    bname="${name%.*}"
    $EQ3NR cmp.d1 $name >/dev/null 2>&1
    for ext in "${EXTENSIONS[@]}"; do
        if ! [ -f $bname.$ext ]; then
            failures=$((failures + 1))
            failedtests+=("$bname.$ext is missing")
        fi
        tests=$((tests + 1))
    done

    rm -rf *

    mkdir local
    cp ../data/cmp/cmp.d1 .
    cp ../data/cmp/acidmwb.3i local/$name
    bname="${name%.*}"
    $EQ3NR cmp.d1 local/$name >/dev/null 2>&1
    for ext in "${EXTENSIONS[@]}"; do
        if ! [ -f $bname.$ext ]; then
            failures=$((failures + 1))
            failedtests+=("$bname.$ext is missing")
        fi
        tests=$((tests + 1))
    done
    cd ..
done

echo "Running tests for $EQ3NR"
for dir in data/*; do
    dataset=$(basename $dir)
    for threei in $dir/*.3i; do
        rm -rf tmp
        mkdir tmp

        problem=$(basename $threei .3i)
        cp $dir/$dataset.d1 tmp
        cp $threei tmp

        cd tmp

        $EQ3NR $dataset.d1 $problem.3i >/dev/null 2>&1
        for file in "${EXTENSIONS[@]}"; do
            file=$problem.$ext
            if [ $file == "$problem.3o" ]; then
                groundtruth=$(realpath ../$dir/$problem.3o)
                assert_same_contents $groundtruth $file skiptimes
            elif [ $file == "$problem.3p" ]; then
                groundtruth=$(realpath ../$dir/$problem.3p)
                assert_same_contents $groundtruth $file
            fi
            tests=$((tests + 1))
        done

        cd ..
        rm -rf tmp
    done
done

summarize
