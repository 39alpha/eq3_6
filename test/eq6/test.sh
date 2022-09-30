#!/usr/bin/env bash

set -e

source $TEST_ROOT/util.sh

EXTENSIONS=(6o 6p)

EQ6=$(realpath -m $PROJECT_ROOT/build/src/eq6/eq6)
if ! [ -f $EQ6 ]; then
    echo "eq6 binary not found at $EQ6; perhaps you should run make?"
    exit -1
fi

echo "  [TEST] Correctly names files"
for name in crisqtz.6i crisqtz; do
    rm -rf tmp
    mkdir tmp
    cd tmp

    cp ../data/cmp/cmp.d1 .
    cp ../data/cmp/crisqtz.6i $name

    bname="${name%.*}"
    $EQ6 cmp.d1 $name >/dev/null
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
    cp ../data/cmp/crisqtz.6i local/$name
    bname="${name%.*}"
    $EQ6 cmp.d1 local/$name >/dev/null
    for ext in "${EXTENSIONS[@]}"; do
        if ! [ -f $bname.$ext ]; then
            failures=$((failures + 1))
            failedtests+=("$bname.$ext is missing")
        fi
        tests=$((tests + 1))
    done
    cd ..
done

echo "  [TEST] Generates correct output"
for dir in data/*; do
    dataset=$(basename $dir)
    for sixi in $dir/*.6i; do
        rm -rf tmp
        mkdir tmp

        problem=$(basename $sixi .6i)
        cp $dir/$dataset.d1 tmp
        cp $sixi tmp

        cd tmp

        $EQ6 $dataset.d1 $problem.6i >/dev/null

        for file in "${EXTENSIONS[@]}"; do
            if [ $file == "$problem.6o" ]; then
                groundtruth=$(realpath ../$dir/$problem.6o)
                assert_same_contents $groundtruth $file skiptimes
            elif [ $file == "$problem.6p" ]; then
                groundtruth=$(realpath ../$dir/$problem.6p)
                assert_same_contents $groundtruth $file
            fi
            tests=$((tests + 1))
        done

        cd ..
        rm -rf tmp
    done
done

summarize
