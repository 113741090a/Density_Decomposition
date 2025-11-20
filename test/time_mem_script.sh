GPP_CMD="g++"
PYTHON_CMD="python3"

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
$PYTHON_CMD plot_single.py
mkdir -p ../../../images/time_mem/mem_normal
mv figures/* ../../../images/time_mem/mem_normal



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
$PYTHON_CMD plot_single.py
mkdir -p ../../../images/time_mem/mem_hyper
mv figures/* ../../../images/time_mem/mem_hyper



cd ../..

chmod +x run_program_time.sh 

./run_program_time.sh

cd normal/data_normal/

$PYTHON_CMD plot_single.py

mkdir -p ../../../images/time_mem/time_normal
mv figures/* ../../../images/time_mem/time_normal


cd ../..

cd hyper/data_hyper/

$PYTHON_CMD plot_single.py


mkdir -p ../../../images/time_mem/time_hyper
mv figures/* ../../../images/time_mem/time_hyper