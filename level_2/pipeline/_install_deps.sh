#!/usr/bin/bash

# download and uncompress julia
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.2-linux-x86_64.tar.gz
tar xvf julia-1.10.2-linux-x86_64.tar.gz

# clone and build TreeQMC
git clone https://github.com/molloy-lab/TREE-QMC.git
cd TREE-QMC/external/MQLib
make
cd ../..
g++ -std=c++11 -O2 \
    -I external/MQLib/include -I external/toms743 \
    -o treeqmc \
    src/*.cpp external/toms743/toms743.cpp \
    external/MQLib/bin/MQLib.a -lm 
