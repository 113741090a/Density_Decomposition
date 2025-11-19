cd parameters/normal/facebook

BOOST_PATH="../../Boost/boost_1_86_0"
#DATASET_URL="https://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz"
WORKING_DIR=$(pwd)
PYTHON_CMD="python3"
GPP_CMD="g++"

# 下载数据集
wget https://snap.stanford.edu/data/facebook_combined.txt.gz -O custom_name.gz
gunzip -c custom_name.gz > graph.txt


# 转换为二分图表示
$PYTHON_CMD norm_to_bi.py


if [ ! -f decomposition ]; then
	$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
fi
./decomposition


mkdir -p data_normal

$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_parameter iterative_normal_parameters.cpp

./iterative_parameter

mkdir -p figures

$PYTHON_CMD draw_all_normal.py

mkdir -p ../../../images/parameters/facebook

mv figures ../../../images/parameters/facebook

