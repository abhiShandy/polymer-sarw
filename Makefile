JDFTX_BUILD_DIR=$(HOME)/Programs/jdftx/build
JDFTx_DIR=$(HOME)/Programs/jdftx/jdftx-git/jdftx

default:
	g++ -std=c++11 -O3 \
	-I $(JDFTx_DIR) -L $(JDFTX_BUILD_DIR) -Wl,-rpath,$(JDFTX_BUILD_DIR) \
	-o main main.cpp -ljdftx

