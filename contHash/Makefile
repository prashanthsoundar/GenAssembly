
ContainmentHash: contHash.o BloomFilter.o MurmurHash3.o MinhashNaive.o
	g++ contHash.o BloomFilter.o MurmurHash3.o MinhashNaive.o -o containmentHash

contHash.o: contHash.cpp
	g++ -c contHash.cpp

BloomFilter.o: BloomFilter.cpp include/BloomFilter.h
	g++ -c BloomFilter.cpp

MinhashNaive.o: MinhashNaive.cpp include/MinhashNaive.h
	g++ -c MinhashNaive.cpp

MurmurHash3.o: MurmurHash3.cpp include/MurmurHash3.h
	g++ -c MurmurHash3.cpp

clean:
	rm *.o containmentHash
