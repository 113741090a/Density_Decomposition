
# Density Decomposition Automation Script
# Version 1.0

# Configuration variables
BOOST_PATH="../../../Boost/boost_1_86_0"
#DATASET_URL="https://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz"
WORKING_DIR=$(pwd)
PYTHON_CMD="python3"
GPP_CMD="g++"

# Define arrays
DATASET_URLS=(
    "https://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz"
    "https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz"
    "https://snap.stanford.edu/data/roadNet-PA.txt.gz"
    "https://snap.stanford.edu/data/roadNet-CA.txt.gz"
    "https://snap.stanford.edu/data/web-Google.txt.gz"
    "https://snap.stanford.edu/data/cit-Patents.txt.gz"
    "https://snap.stanford.edu/data/wiki-topcats.txt.gz"
    "https://snap.stanford.edu/data/bigdata/communities/com-orkut.ungraph.txt.gz"
    # You can add more datasets here. 
)

# The original data downloaded from the internet may contain unnecessary lines. 
# Here each number corresponds to a dataset, meaning the number of unnecessary lines that need to be deleted. 
PYTHON_PARAMS=(
    "4"  
    "0"
    "5"
    "0"
    "0"
    "4"
    "0"
    "4"
)

# The number of iterations for normal graph
ITERATION_COUNTS1=(
    "1400"  
    "750"
    "3000"
    "4200"
    "1800"
    "4100"
    "2800"
    "1000"
)

# The number of iterations for hypergraph
ITERATION_COUNTS2=(
    "400"  # Iteration count for the first dataset
    "750"
    "3000"
    "700"
    "1200"
    "2000"
    "1500"
    "200"
    # You can add more iteration counts
)

# Folder name where the data for each dataset is stored. 
OUTPUT_PATHS=(
    "amazon"  
    "dblp"
    "Roadnet-PA"
    "Roadnet-CA"
    "web-google"
    "patents"
    "wiki-cats"
    "orkut"
    # You can add more output paths
)

BIPARTITE_NAMES=(
    "amazon.txt"  # New name for the first dataset's bipartite_hyper.txt
    "dblp.txt"
    "roadnet-PA.txt"
    "roadnet-CA.txt"
    "web-google.txt"
    "patent.txt"
    "top-cats.txt"
    "orkut.txt"
)

HYPER_NAMES=(
    "amazon.txt"  # New name for the first dataset's bipartite_hyper.txt
    "dblp.txt"
    "roadnet-PA.txt"
    "roadnet-CA.txt"
    "web-google.txt"
    "patent.txt"
    "top-cats.txt"
    "orkut.txt"
    # You can add more names
)

WORKING_DIR=$(pwd)
PYTHON_CMD="python3"
GPP_CMD="g++"

# Check if array lengths are consistent
if [ ${#DATASET_URLS[@]} -ne ${#PYTHON_PARAMS[@]} ] || 
   [ ${#DATASET_URLS[@]} -ne ${#ITERATION_COUNTS1[@]} ] || 
   [ ${#DATASET_URLS[@]} -ne ${#OUTPUT_PATHS[@]} ] ||
   [ ${#DATASET_URLS[@]} -ne ${#BIPARTITE_NAMES[@]} ] ||
   [ ${#DATASET_URLS[@]} -ne ${#HYPER_NAMES[@]} ]; then
    echo "Error: Array lengths do not match!"
    echo "Number of dataset URLs: ${#DATASET_URLS[@]}"
    echo "Number of Python parameters: ${#PYTHON_PARAMS[@]}"
    echo "Number of iteration counts: ${#ITERATION_COUNTS1[@]}"
    echo "Number of output paths: ${#OUTPUT_PATHS[@]}"
    echo "Number of bipartite file names: ${#BIPARTITE_NAMES[@]}"
    echo "Number of hypergraph file names: ${#HYPER_NAMES[@]}"
    exit 1
fi


# Create directory to save processed results
mkdir -p processed_graphs

# Process all datasets
for i in "${!DATASET_URLS[@]}"; do
    DATASET_URL="${DATASET_URLS[$i]}"
    PYTHON_PARAM="${PYTHON_PARAMS[$i]}"
    ITERATION_COUNT1="${ITERATION_COUNTS1[$i]}"
    ITERATION_COUNT2="${ITERATION_COUNTS2[$i]}"
    OUTPUT_PATH="${OUTPUT_PATHS[$i]}"
    BIPARTITE_NAME="${BIPARTITE_NAMES[$i]}"
    HYPER_NAME="${HYPER_NAMES[$i]}"
    
    echo "Processing dataset $((i+1))/${#DATASET_URLS[@]}: $DATASET_URL"
    echo "Parameters: $PYTHON_PARAM"
    echo "Iteration count: $ITERATION_COUNT1"
    echo "Iteration count: $ITERATION_COUNT2"
    echo "Output path: $OUTPUT_PATH"
    echo "Bipartite file will be saved as: $BIPARTITE_NAME"
    echo "Hypergraph file will be saved as: $HYPER_NAME"
    
    # Data preprocessing
    echo "Starting data preprocessing..."
    cd Density_Decomposition
    cd data-preprocessing

    # Download dataset
    wget $DATASET_URL -O custom_name.gz
    gunzip -c custom_name.gz > graph.txt

    # Process graph file
    $PYTHON_CMD process.py $PYTHON_PARAM

    # Convert to bipartite representation
    $PYTHON_CMD norm_to_bi.py
    $PYTHON_CMD norm_to_hyper.py

    echo "Saving processed graph files..."
    cp bipartite_mark.txt "../time_mem/normal/${BIPARTITE_NAME}"
    cp bipartite_hyper.txt "../time_mem/hyper/${HYPER_NAME}"

    # Move processed files to corresponding directories
    mv bipartite_mark.txt ../normal/
    mv bipartite_hyper.txt ../hyper/

    cd ..

    # Clean temporary files
    rm -f custom_name.gz
    
    #cd ..

    # Exact decomposition - Normal graph
    echo "Running exact decomposition for normal graph..."
    cd normal
    # Compile only once (can skip if already compiled)
    if [ ! -f decomposition ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
    fi
    ./decomposition
    cd ..

    # Exact decomposition - Hypergraph
    echo "Running exact decomposition for hypergraph..."
    cd hyper
    # Compile only once (can skip if already compiled)
    if [ ! -f decomposition ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
    fi
    ./decomposition
    cd ..

    # Approximation algorithm - Normal graph
    echo "Running approximation algorithm for normal graph..."
    cd normal
    mkdir -p data_normal
    # Compile only once (can skip if already compiled)
    if [ ! -f iterative_normal ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_normal iterative_normal.cpp
    fi
    ./iterative_normal $ITERATION_COUNT1
    cd ..

    # Approximation algorithm - Hypergraph
    echo "Running approximation algorithm for hypergraph..."
    cd hyper
    mkdir -p data_hyper
    # Compile only once (can skip if already compiled)
    if [ ! -f iterative_hyper ]; then
        $GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_hyper iterative_hyper.cpp
    fi
    ./iterative_hyper $ITERATION_COUNT2
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
    mkdir -p ../images/$OUTPUT_PATH
    mv normal/figures ../images/$OUTPUT_PATH/figures_normal
    mv hyper/figures ../images/$OUTPUT_PATH/figures_hyper
    
    echo "Completed processing dataset $((i+1))/${#DATASET_URLS[@]}"
    echo "----------------------------------------"
done

echo "All datasets processed successfully!"
echo "Processed graph files are saved in: ${WORKING_DIR}/processed_graphs/"

# Time and Memory Measurements - Normal Graph
echo "Running time and memory tests for normal graph..."
cd time_mem/normal
mkdir -p data_normal data_normal_mem figures
$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o test_time test_time.cpp
$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o test_mem test_memory.cpp

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