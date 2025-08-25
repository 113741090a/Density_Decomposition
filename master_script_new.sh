#!/bin/bash

# Density Decomposition Automation Script
# Version 1.0

# Configuration variables
BOOST_PATH="../../Boost/boost_1_86_0"
#DATASET_URL="https://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz"
WORKING_DIR=$(pwd)
PYTHON_CMD="python3"
GPP_CMD="g++"

# Create directory structure
# mkdir -p data-preprocessing \
#          normal \
#          hyper \
#          time_mem/normal \
#          time_mem/hyper \
#          figures

#unzip Boost/boost_1_86_0
cat Boost.zip.part* > Boost.zip




# Install dependencies
echo "Installing dependencies..."
# sudo apt-get update
# sudo apt-get install -y python3-pip g++ git wget unzip
pip install matplotlib




# 定义数组
DATASET_URLS=(
    "https://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz"
    "https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz"
    "https://snap.stanford.edu/data/roadNet-PA.txt.gz"
    "https://snap.stanford.edu/data/roadNet-CA.txt.gz"
    "https://snap.stanford.edu/data/web-Google.txt.gz"
    "https://snap.stanford.edu/data/cit-Patents.txt.gz"
    "https://snap.stanford.edu/data/wiki-topcats.txt.gz"
    "https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz"
    # 可以添加更多数据集URL
)

PYTHON_PARAMS=(
    "4"  # 对应第一个数据集的参数
    "0"
    "5"
    "0"
    "0"
    "4"
    "0"
    "4"
    # 可以添加更多参数
)

ITERATION_COUNTS1=(
    "1400"  # 第一个数据集的迭代次数
    "750"
    "3000"
    "4200"
    "1800"
    "4100"
    "2800"
    "1000"
    # 可以添加更多迭代次数
)

ITERATION_COUNTS2=(
    "400"  # 第一个数据集的迭代次数
    "750"
    "3000"
    "700"
    "1200"
    "2000"
    "1500"
    "200"
    # 可以添加更多迭代次数
)

OUTPUT_PATHS=(
    "amazon"  # 第一个数据集的输出路径
    "dblp"
    "Roadnet-PA"
    "Roadnet-CA"
    "web-google"
    "patents"
    "wiki-cats"
    "orkut"
    # 可以添加更多输出路径
)

BIPARTITE_NAMES=(
    "amazon.txt"  # 第一个数据集bipartite_hyper.txt的新名称
    "dblp.txt"
    "roadnet-PA.txt"
    "roadnet-CA.txt"
    "web-google.txt"
    "patent.txt"
    "top-cats.txt"
    "orkut.txt"
)

HYPER_NAMES=(
    "amazon.txt"  # 第一个数据集bipartite_hyper.txt的新名称
    "dblp.txt"
    "roadnet-PA.txt"
    "roadnet-CA.txt"
    "web-google.txt"
    "patent.txt"
    "top-cats.txt"
    "orkut.txt"
    #("amazon.txt" "dblp.txt" "roadnet-PA.txt" "roadnet-CA.txt" "web-google.txt" "patent.txt" "top-cats.txt" "orkut.txt")
    # 可以添加更多新名称
)

WORKING_DIR=$(pwd)
PYTHON_CMD="python3"
GPP_CMD="g++"

# 检查数组长度是否一致
if [ ${#DATASET_URLS[@]} -ne ${#PYTHON_PARAMS[@]} ] || 
   [ ${#DATASET_URLS[@]} -ne ${#ITERATION_COUNTS1[@]} ] || 
   [ ${#DATASET_URLS[@]} -ne ${#OUTPUT_PATHS[@]} ] ||
   [ ${#DATASET_URLS[@]} -ne ${#BIPARTITE_NAMES[@]} ] ||
   [ ${#DATASET_URLS[@]} -ne ${#HYPER_NAMES[@]} ]; then
    echo "错误：数组长度不一致！"
    echo "数据集URL数量: ${#DATASET_URLS[@]}"
    echo "Python参数数量: ${#PYTHON_PARAMS[@]}"
    echo "迭代次数数量: ${#ITERATION_COUNTS1[@]}"
    echo "输出路径数量: ${#OUTPUT_PATHS[@]}"
    echo "二分图文件名称数量: ${#BIPARTITE_NAMES[@]}"
    echo "超图文件名称数量: ${#HYPER_NAMES[@]}"
    exit 1
fi

# 创建保存处理结果的目录
mkdir -p processed_graphs

# 遍历所有数据集
for i in "${!DATASET_URLS[@]}"; do
    DATASET_URL="${DATASET_URLS[$i]}"
    PYTHON_PARAM="${PYTHON_PARAMS[$i]}"
    ITERATION_COUNT1="${ITERATION_COUNTS1[$i]}"
    ITERATION_COUNT2="${ITERATION_COUNTS2[$i]}"
    OUTPUT_PATH="${OUTPUT_PATHS[$i]}"
    BIPARTITE_NAME="${BIPARTITE_NAMES[$i]}"
    HYPER_NAME="${HYPER_NAMES[$i]}"
    
    echo "处理数据集 $((i+1))/${#DATASET_URLS[@]}: $DATASET_URL"
    echo "参数: $PYTHON_PARAM"
    echo "迭代次数: $ITERATION_COUNT1"
    echo "迭代次数: $ITERATION_COUNT2"
    echo "输出路径: $OUTPUT_PATH"
    echo "二分图文件将保存为: $BIPARTITE_NAME"
    echo "超图文件将保存为: $HYPER_NAME"
    
    # 数据预处理
    echo "开始数据预处理..."
    cd Density_Decomposition
    cd data-preprocessing

    # 下载数据集
    wget $DATASET_URL -O custom_name.gz
    gunzip -c custom_name.gz > graph.txt

    # 处理图文件
    $PYTHON_CMD process.py $PYTHON_PARAM

    # 转换为二分图表示
    $PYTHON_CMD norm_to_bi.py
    $PYTHON_CMD norm_to_hyper.py

    echo "保存处理后的图文件..."
    cp bipartite_mark.txt "../time_mem/normal/${BIPARTITE_NAME}"
    cp bipartite_hyper.txt "../time_mem/hyper/${HYPER_NAME}"

    # 移动处理后的文件到相应目录
    mv bipartite_mark.txt ../normal/
    mv bipartite_hyper.txt ../hyper/
    # 保存处理后的文件到指定位置


    cd ..



    
    
    # 清理临时文件
    rm -f custom_name.gz
    
    cd ..

    # 精确分解 - 普通图
    echo "运行普通图的精确分解..."
    cd normal
    # 只编译一次（如果已经编译过可以跳过）
    if [ ! -f decomposition ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
    fi
    ./decomposition
    cd ..

    # 精确分解 - 超图
    echo "运行超图的精确分解..."
    cd hyper
    # 只编译一次（如果已经编译过可以跳过）
    if [ ! -f decomposition ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
    fi
    ./decomposition
    cd ..

    # 近似算法 - 普通图
    echo "运行普通图的近似算法..."
    cd normal
    mkdir -p data_normal
    # 只编译一次（如果已经编译过可以跳过）
    if [ ! -f iterative_normal ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_normal iterative_normal.cpp
    fi
    ./iterative_normal $ITERATION_COUNT1
    cd ..

    # 近似算法 - 超图
    echo "运行超图的近似算法..."
    cd hyper
    mkdir -p data_hyper
    # 只编译一次（如果已经编译过可以跳过）
    if [ ! -f iterative_hyper ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_hyper iterative_hyper.cpp
    fi
    ./iterative_hyper $ITERATION_COUNT2
    cd ..

    # 数据可视化 - 普通图
    echo "生成普通图的可视化..."
    cd normal
    mkdir -p figures
    $PYTHON_CMD draw_all_normal.py
    cd ..

    # 数据可视化 - 超图
    echo "生成超图的可视化..."
    cd hyper
    mkdir -p figures
    $PYTHON_CMD draw_all_hyper.py
    cd ..

    # 创建输出目录并移动结果
    mkdir -p ../images/$OUTPUT_PATH
    mv normal/figures ../images/$OUTPUT_PATH/figures_normal
    mv hyper/figures ../images/$OUTPUT_PATH/figures_hyper
    
    echo "完成数据集 $((i+1))/${#DATASET_URLS[@]} 的处理"
    echo "----------------------------------------"
done

echo "所有数据集处理完成！"
echo "处理后的图文件已保存在: ${WORKING_DIR}/processed_graphs/"



# Time and Memory Measurements - Normal Graph
echo "Running time and memory tests for normal graph..."
cd time_mem/normal
mkdir -p data_normal data_normal_mem figures
$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o test_time test_time.cpp
$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o test_mem test_memory.cpp

# Update params in run_program_time.sh and run_memory.sh as needed
# Here we use a default set of parameters
# echo 'params=("bipartite_mark.txt")' > run_program_time.sh
# echo 'params=("bipartite_mark.txt")' > run_memory.sh
# echo 'algos=("Elist++" "FISTA*" "Greedy++" "PR" "PR_EXP" "PR_LIN")' >> run_memory.sh

chmod +x run_memory.sh
./run_memory.sh

cd data_normal_mem/
$PYTHON_CMD plot.py
mkdir -p ../../../../images/time_mem
mv figures ../../../../images/time_mem/mem_normal

cd ../

cd ../..


# Time and Memory Measurements - Hypergraph
echo "Running time and memory tests for hypergraph..."
cd time_mem/hyper
mkdir -p data_hyper data_hyper_mem figures
$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o test_time_hyper test_time_hyper.cpp
$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o test_mem_hyper test_mem_hyper.cpp

# Run test memory
chmod +x run_memory.sh
./run_memory.sh



cd data_hyper_mem/
$PYTHON_CMD plot.py
#mkdir -p ../../../images/time_mem
mv figures ../../../../images/time_mem/mem_hyper



cd ../..

chmod +x run_program_time.sh

./run_program_time.sh

cd normal/data_normal/

$PYTHON_CMD plot.py

mv figures ../../../../images/time_mem/time_normal


cd ../..

cd hyper/data_hyper/

$PYTHON_CMD plot.py

mv figures ../../../../images/time_mem/time_hyper

cd ../..

echo "All tasks completed successfully!"
