default:
	g++ -std=c++11 -O3 \
	-I $(JDFTx_DIR) \
	-o main main.cpp
