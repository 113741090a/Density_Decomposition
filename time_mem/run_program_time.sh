#!/bin/bash


# datasets names
params=("patent.txt" "top-cats.txt" "orkut.txt")
# locations of your program
normal_program="./normal/test_time"
hyper_program="./hyper/test_time_hyper"


# execute 
echo "Running normal program:"
for param in "${params[@]}"; do
    echo "Executing $normal_program with parameter: $param"
    $normal_program "normal/$param"
done


# 执行 hyper 程序
echo "Running hyper program:"
for param in "${params[@]}"; do
    echo "Executing $hyper_program with parameter: $param"
    $hyper_program "hyper/$param"
done

