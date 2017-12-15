#include "include/BloomFilter.h"
#include<iostream>
#include<math.h>
#include "include/MurmurHash3.h"
#include<stdlib.h>

using namespace std;

BloomFilter::BloomFilter(uint64_t capacity,double errorRate)
{
    this->capacity=capacity;
    this->errorRate=errorRate;
    BloomFilter::initBloom(capacity,errorRate);
}

void BloomFilter::printSetBits()
{
    for(int i=0;i<10000000;i++)
    {
        if(m_bits[i])
            cout<<i<<" ";
    }
}
void BloomFilter::initBloom(uint64_t capacity,double errorRate)
{
    uint64_t filterSize;
    uint64_t numHashes;
    
    double num = capacity*abs(log(errorRate));
    double denom = log(2)*log(2);
    
    filterSize = (uint64_t)ceil(num/denom);
   // this->m_bits.resize(filterSize);
    
    double part1 = filterSize/capacity;
    numHashes = (uint64_t) ceil(part1*log(2));
    this->numHashes = numHashes;
}

uint64_t BloomFilter::hash(string data,uint64_t filterSize,int n) const{
    uint64_t hashVal[2];
    MurmurHash3_x64_128(data.c_str(), data.length(), 0, &hashVal);
    return (hashVal[0] + n * hashVal[1]) % filterSize;
}

void BloomFilter::add(string data){
    
    for (int n = 1; n <= numHashes; n++) {
        m_bits[hash(data,10000000,n)] = true;
    }
}

uint64_t BloomFilter::getMBitsSize()
{
    return 10000000;
}

uint64_t BloomFilter::getNumHashes()
{
    return numHashes;
}
uint64_t BloomFilter::getBSize(){
    return b_size;
}

void BloomFilter::setBSize(uint64_t size){
    this->b_size  = size;
}
bool BloomFilter::possiblyContains(string data) const {
    for (int n = 1; n <= numHashes; n++) {
        if (!m_bits[hash(data,10000000,n)]) {
            return false;
        }
    }
    
    return true;
}

