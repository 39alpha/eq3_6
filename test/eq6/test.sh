#!/usr/bin/env bash

source $TEST_ROOT/util.sh

EQ6=$(realpath -m $PROJECT_ROOT/bin/eq6)
if ! [ -f $EQ6 ]; then
    echo "eq6 binary not found at $EQ6; perhaps you should run make?"
    exit -1
fi

echo "Running tests for $EQ6"
for dir in data/*; do
    dataset=$(basename $dir)
    for sixi in $dir/*.6i; do
        rm -rf tmp
        mkdir tmp

        problem=$(basename $sixi .6i)
        cp $dir/data1.$dataset tmp
        cp $sixi tmp

        cd tmp

        $EQ6 data1.$dataset $problem.6i >/dev/null 2>&1

        for file in output pickup; do
            if [ $file == "output" ]; then
                groundtruth=$(realpath ../$dir/$problem.6o)
                assert_same_contents $groundtruth $file skiptimes
            elif [ $file == "pickup" ]; then
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
