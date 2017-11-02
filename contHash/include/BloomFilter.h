#ifndef _bloom_h
#define _bloom_h
#include<vector>
#include<string>
using namespace std;

typedef unsigned long long int uint64_t;

class BloomFilter
{
    uint64_t capacity;
    double errorRate;

    int numHashes;
    vector<bool> m_bits;

    public:
        BloomFilter(uint64_t capacity,double errorRate);
        void initBloom(uint64_t capacity,double errorRate);
        void add(string data);
        void printSetBits();
        bool possiblyContains(string data) const;
        uint64_t hash(string data,uint64_t filterSize,int n) const;
        ~BloomFilter();

};

#endif
