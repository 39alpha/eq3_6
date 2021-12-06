#!/usr/bin/env bash

set -e

source $TEST_ROOT/util.sh

EXTENSIONS=(po d1 d1f s)

EQPT=$(realpath -m $PROJECT_ROOT/bin/eqpt)
if ! [ -f $EQPT ]; then
    echo "eqpt binary not found at $EQPT; perhaps you should run make?"
    exit -1
fi

echo "Running tests for $EQPT"

echo "  [TEST] Correctly names files"
for name in data0.1kb 1kb.d0 1kb; do
    rm -rf tmp
    mkdir tmp
    cd tmp

    cp ../data/1kb.d0 $name
    bname="${name%.*}"
    $EQPT $name >/dev/null
    for ext in "${EXTENSIONS[@]}"; do
        if ! [ -f $bname.$ext ]; then
            failures=$((failures + 1))
            failedtests+=("$bname.$ext is missing")
        fi
        tests=$((tests + 1))
    done

    rm -rf *

    mkdir local
    cp ../data/1kb.d0 local/$name
    bname="${name%.*}"
    $EQPT local/$name >/dev/null
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
for data0 in data/*.d0; do
    dataset=$(basename $data0 .d0)

    rm -rf tmp
    mkdir tmp

    cp $data0 tmp

    cd tmp
    $EQPT $(basename $data0) >/dev/null
    for ext in "${EXTENSIONS[@]}"; do
        file=$dataset.$ext
        groundtruth=$(realpath ../data/$dataset.$ext)
        if [ $file == "$dataset.po" ]; then
            assert_same_contents $groundtruth $file skiptimes
        else
            assert_same_contents $groundtruth $file
        fi
        tests=$((tests + 1))
    done

    cd ..
    rm -rf tmp
done

summarize
