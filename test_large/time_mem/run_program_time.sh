#!/bin/bash


# datasets names
#params=("amazon.txt" "dblp.txt" "roadnet-PA.txt" "roadnet-CA.txt" "web-google.txt" "patent.txt" "top-cats.txt" "orkut.txt") # "orkut.txt"

params=("top-cats.txt") # "amazon.txt" "web-google.txt" "orkut.txt"  "roadnet-PA.txt" "roadnet-CA.txt" "web-google.txt" "patent.txt" "top-cats.txt" "orkut.txt"
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

