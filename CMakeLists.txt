cmake_minimum_required(VERSION 3.9)
project(GenAssembly)

set(CMAKE_CXX_STANDARD 11)

add_executable(GenAssembly
        contHash/include/BloomFilter.h
        contHash/include/MinhashNaive.h
        contHash/include/MurmurHash3.h
        contHash/BloomFilter.cpp
        contHash/contHash.cpp
        contHash/MinhashNaive.cpp
        contHash/MurmurHash3.cpp
        contHash/README.md
        test/GenAssembly/contHash/include/BloomFilter.h
        test/GenAssembly/contHash/include/MinhashNaive.h
        test/GenAssembly/contHash/include/MurmurHash3.h
        test/GenAssembly/contHash/BloomFilter.cpp
        test/GenAssembly/contHash/contHash.cpp
        test/GenAssembly/contHash/MinhashNaive.cpp
        test/GenAssembly/contHash/MurmurHash3.cpp)
