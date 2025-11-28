#!/usr/bin/env bash

# Configuration variables
BOOST_PATH="../Boost/boost_1_86_0"

PYTHON_CMD="python3"

GPP_CMD="g++"




# Folder name where the data for each dataset is stored. 
OUTPUT_PATHS=(
    "facebook"
    # You can add more output paths
)


BIPARTITE_NAMES=(
    "facebook.txt"
)


HYPER_NAMES=(
    "facebook.txt"
)

# The number of iterations for normal graph
ITERATION_COUNTS1=(
    "800"
)

# The number of iterations for hypergraph
ITERATION_COUNTS2=(
    "800"
    # You can add more iteration counts
)



# Process all datasets
for i in "${!OUTPUT_PATHS[@]}"; do

    
    #Set names
    ITERATION_COUNT1="${ITERATION_COUNTS1[$i]}"
    ITERATION_COUNT2="${ITERATION_COUNTS2[$i]}"
    OUTPUT_PATH="${OUTPUT_PATHS[$i]}"
    BIPARTITE_NAME="${BIPARTITE_NAMES[$i]}"
    HYPER_NAME="${HYPER_NAMES[$i]}"

    # Exact decomposition - Normal graph
    echo "Running exact decomposition for normal graph for graph ..."
    cd normal
    # Compile only once (can skip if already compiled)
    if [ ! -f decomposition ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
    fi
    ./decomposition $BIPARTITE_NAME
    cd ..

    # Exact decomposition - Hypergraph
    echo "Running exact decomposition for hypergraph..."
    cd hyper
    # Compile only once (can skip if already compiled)
    if [ ! -f decomposition ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
    fi
    ./decomposition $HYPER_NAME
    cd ..

    # Approximation algorithm - Normal graph
    echo "Running approximation algorithm for normal graph..."
    cd normal
    dir_path="data_normal"
    if [ -d "$dir_path" ]; then
        echo "即将删除: $dir_path"
        rm -r "$dir_path"
    fi
    mkdir -p data_normal
    # Compile only once (can skip if already compiled)
    if [ ! -f iterative_normal ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_normal iterative_normal.cpp
    fi
    ./iterative_normal $BIPARTITE_NAME $ITERATION_COUNT1
    cd ..

    # Approximation algorithm - Hypergraph
    echo "Running approximation algorithm for hypergraph..."
    cd hyper
    dir_path="data_hyper"
    if [ -d "$dir_path" ]; then
        echo "即将删除: $dir_path"
        rm -r "$dir_path"
    fi
    mkdir -p data_hyper
    # Compile only once (can skip if already compiled)
    if [ ! -f iterative_hyper ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_hyper iterative_hyper.cpp
    fi
    ./iterative_hyper $HYPER_NAME $ITERATION_COUNT2
    cd ..

    # Data visualization - Normal graph
    echo "Generating visualization for normal graph..."
    cd normal
    mkdir -p figures
    $PYTHON_CMD draw_all_normal.py
    cd ..

    # Data visualization - Hypergraph
    echo "Generating visualization for hypergraph..."
    cd hyper
    mkdir -p figures
    $PYTHON_CMD draw_all_hyper.py
    cd ..

    # Create output directory and move results
    #images地址换了。
    mkdir -p images/$OUTPUT_PATH
    mkdir -p images/$OUTPUT_PATH/figures_normal
    mkdir -p images/$OUTPUT_PATH/figures_hyper

    mv normal/figures/* images/$OUTPUT_PATH/figures_normal
    mv hyper/figures/* images/$OUTPUT_PATH/figures_hyper

    # rm -r normal/data_normal

    # rm -r hyper/data_hyper


    
    echo "Completed processing dataset $((i+1))/${#DATASET_URLS[@]}"
    echo "----------------------------------------"



done
