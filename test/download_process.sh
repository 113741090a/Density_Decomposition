#python command name
PYTHON_CMD="python3"


# Define arrays
#https://snap.stanford.edu/data/facebook_combined.txt.gz
DATASET_URLS=(
    "https://snap.stanford.edu/data/facebook_combined.txt.gz"
    # You can add more datasets here. 
)

# The original data downloaded from the internet may contain unnecessary lines. 
# Here each number corresponds to a dataset, meaning the number of unnecessary lines that need to be deleted. 
PYTHON_PARAMS=(
    "0"
)


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

# Check if array lengths are consistent
if [ ${#DATASET_URLS[@]} -ne ${#PYTHON_PARAMS[@]} ] || 
   [ ${#DATASET_URLS[@]} -ne ${#OUTPUT_PATHS[@]} ] ||
   [ ${#DATASET_URLS[@]} -ne ${#BIPARTITE_NAMES[@]} ] ||
   [ ${#DATASET_URLS[@]} -ne ${#HYPER_NAMES[@]} ]; then
    echo "Error: Array lengths do not match!"
    echo "Number of dataset URLs: ${#DATASET_URLS[@]}"
    echo "Number of Python parameters: ${#PYTHON_PARAMS[@]}"
    echo "Number of output paths: ${#OUTPUT_PATHS[@]}"
    echo "Number of bipartite file names: ${#BIPARTITE_NAMES[@]}"
    echo "Number of hypergraph file names: ${#HYPER_NAMES[@]}"
    exit 1
fi

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
    #cd Density_Decomposition
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

    #This is used for later testing for time and memory
    cp bipartite_mark.txt "../time_mem/normal/${BIPARTITE_NAME}"
    cp bipartite_hyper.txt "../time_mem/hyper/${HYPER_NAME}"

    cp bipartite_mark.txt "../normal/${BIPARTITE_NAME}"
    cp bipartite_hyper.txt "../hyper/${HYPER_NAME}"

    # Clean temporary files
    rm -f custom_name.gz
    rm graph.txt
    rm bipartite_mark.txt
    rm bipartite_hyper.txt


done