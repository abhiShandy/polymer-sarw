JDFTX_BUILD_DIR=$(HOME)/Code/JDFTx/build

default:
	g++ -std=c++11 -O3 \
	-I $(JDFTx_DIR) -L $(JDFTX_BUILD_DIR) -Wl,-rpath,$(JDFTX_BUILD_DIR) \
	-o main main.cpp -ljdftx

