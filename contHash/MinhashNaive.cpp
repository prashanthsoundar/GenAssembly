#include "include/MinhashNaive.h"
#include "include/MurmurHash3.h"
#include<climits>
#include<iostream>
using namespace std;

MinhashNaive::MinhashNaive(uint64_t numHashes,uint64_t largeMod)
{
    this->numHashes = numHashes;
    this->largeMod = largeMod;
    for(int i=0;i<numHashes;i++)
    {
        sketch.push_back(ULLONG_MAX);
        sketchStr.push_back(" ");
    }
    
}

uint64_t MinhashNaive::hash(string data,int n) const{
    uint64_t hashVal[2];
    MurmurHash3_x64_128(data.c_str(), data.length(), 0, &hashVal);
    return (hashVal[0] + n * hashVal[1]) % largeMod;
}

vector<string> MinhashNaive::computeHashedKmers(string s,uint64_t kmerSize)
{
    for (int i = 0; i < s.length() - kmerSize+1; i++) {
        string sub_string = s.substr(i, kmerSize);
        for (int j = 0; j < numHashes; j++) {
            uint64_t temp = MinhashNaive::hash(sub_string, j);
            if (sketch[j] > temp) {
                sketch[j] = temp;
                sketchStr[j] = sub_string;
            }
        }
        
    }
    return sketchStr;
    
}
