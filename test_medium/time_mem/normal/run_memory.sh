#!/bin/bash

# 定义要传递的参数列表
#params=("amazon.txt" "dblp.txt" "roadnet-PA.txt" "roadnet-CA.txt" "web-google.txt" "patent.txt" "top-cats.txt" "orkut.txt") #("amazon.txt" "patent.txt" "top-cats.txt")
params=("amazon.txt") #"amazon.txt" "web-google.txt"

algos=("FISTA*" "PR_exp" "Elist++" "PR" "PR_lin" "Greedy++")
#params=("bipartite_mark2.txt" "bipartite_mark3.txt" "bipartite_mark4.txt" "bipartite_mark5.txt")
# 定义程序路径
normal_program="./test_mem"
hyper_program="./test_mem_hyper"

# 执行 normal 程序
echo "Running normal program:"
for param in "${params[@]}"; do
    for algo in "${algos[@]}"; do
        echo "Executing $normal_program with parameter: $param"
        $normal_program "$param" "$algo"
    done
done

# # 执行 hyper 程序
# echo "Running hyper program:"
# for param in "${params[@]}"; do
#     echo "Executing $hyper_program with parameter: $param"
#     $hyper_program "hyper/$param"
# done