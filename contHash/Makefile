
ContainmentHash: contHash.o BloomFilter.o MurmurHash3.o MinhashNaive.o
	g++ -std=c++0x contHash.o BloomFilter.o MurmurHash3.o MinhashNaive.o -o containmentHash

contHash.o: contHash.cpp
	g++ -std=c++0x -c contHash.cpp

BloomFilter.o: BloomFilter.cpp include/BloomFilter.h
	g++  -std=c++0x -c BloomFilter.cpp

MinhashNaive.o: MinhashNaive.cpp include/MinhashNaive.h
	g++  -std=c++0x -c MinhashNaive.cpp

MurmurHash3.o: MurmurHash3.cpp include/MurmurHash3.h
	g++ -std=c++0x -c MurmurHash3.cpp

clean:
	rm *.o containmentHash
