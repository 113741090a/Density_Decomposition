cd parameters/normal/web-google

BOOST_PATH="../../Boost/boost_1_86_0"
#DATASET_URL="https://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz"
WORKING_DIR=$(pwd)
PYTHON_CMD="python3"
GPP_CMD="g++"

# 下载数据集
wget https://snap.stanford.edu/data/web-Google.txt.gz -O custom_name.gz
gunzip -c custom_name.gz > graph.txt


# 转换为二分图表示
$PYTHON_CMD norm_to_bi.py


if [ ! -f decomposition ]; then
	$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
fi
./decomposition


$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_parameter iterative_normal_parameters.cpp

./iterative_parameter

$PYTHON_CMD draw_all_normal.py

mkdir -p ../../../images/parameters/web-google

mv figures ../../../images/parameters/web-google

cd ../
cd wiki-cats/

# 下载数据集
wget https://snap.stanford.edu/data/wiki-topcats.txt.gz -O custom_name.gz
gunzip -c custom_name.gz > graph.txt


# 转换为二分图表示
$PYTHON_CMD norm_to_bi.py

if [ ! -f decomposition ]; then
	$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition decomposition.cpp
fi
./decomposition


$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_parameter iterative_normal_parameters.cpp

./iterative_parameter

$PYTHON_CMD draw_all_normal.py

mkdir -p ../../../images/parameters/wiki-cats

mv figures ../../../images/parameters/wiki-cats


cd ../../../

cd weighted/

BOOST_PATH="../Boost/boost_1_86_0"
$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o decomposition_new decomposition_new.cpp

$GPP_CMD -I$BOOST_PATH -L$BOOST_PATH/libs -o iterative_normal iterative_normal_new.cpp

./iterative_normal


$PYTHON_CMD draw_all_normal.py

mkdir -p ../images/wiki-selec-weighted/figures_normal

mv figures ../images/wiki-selec-weighted/figures_normal

echo "successful!"