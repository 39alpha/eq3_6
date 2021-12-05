#!/usr/bin/env bash

source $TEST_ROOT/util.sh

EQPT=$(realpath -m $PROJECT_ROOT/bin/eqpt)
if ! [ -f $EQPT ]; then
    echo "eqpt binary not found at $EQPT; perhaps you should run make?"
    exit -1
fi

echo "Running tests for $EQPT"
for dir in data/*; do
    dataset=$(basename $dir)
    data0=$dir/data0.$dataset

    rm -rf tmp
    mkdir tmp

    problem=$(basename $data0)
    cp $data0 tmp

    cd tmp
    $EQPT data0.$dataset >/dev/null 2>&1
    for file in output data1 data1f slist; do
        groundtruth=$(realpath ../$dir/$file.$dataset)
        if [ $file == "output" ]; then
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
