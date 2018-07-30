# Add two variables to the environment pointing towards JDFTx source and build directory
# Example:
# export JDFTx_SRC_DIR="$HOME/Code/JDFTx/jdftx-git/jdftx/"
# export JDFTx_BUILD_DIR="$HOME/Code/JDFTx/build/"

default:
	g++ -std=c++11 -O3 \
	-I $(JDFTx_SRC_DIR) -L $(JDFTx_BUILD_DIR) -Wl,-rpath,$(JDFTx_BUILD_DIR) \
	-o main main.cpp -ljdftx
