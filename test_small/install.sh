#!/bin/bash




#This script is to install necessary 

#unzip Boost/boost_1_86_0
cat ../Boost.zip.part* > Boost.zip
unzip Boost.zip

# Install dependencies
echo "Installing dependencies..."
# sudo apt-get update
# sudo apt-get install -y python3-pip g++ git wget unzip
pip install matplotlib