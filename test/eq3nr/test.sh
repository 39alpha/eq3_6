#!/usr/bin/env bash

source $TEST_ROOT/util.sh

EQ3NR=$(realpath -m $PROJECT_ROOT/bin/eq3nr)
if ! [ -f $EQ3NR ]; then
    echo "eq3nr binary not found at $EQ3NR; perhaps you should run make?"
    exit -1
fi

echo "Running tests for $EQ3NR"
for dir in data/*; do
    dataset=$(basename $dir)
    for threei in $dir/*.3i; do
        rm -rf tmp
        mkdir tmp

        problem=$(basename $threei .3i)
        cp $dir/data1.$dataset tmp
        cp $threei tmp

        cd tmp

        $EQ3NR data1.$dataset $problem.3i >/dev/null 2>&1
        for file in output pickup; do
            if [ $file == "output" ]; then
                groundtruth=$(realpath ../$dir/$problem.3o)
                assert_same_contents $groundtruth $file skiptimes
            elif [ $file == "pickup" ]; then
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
