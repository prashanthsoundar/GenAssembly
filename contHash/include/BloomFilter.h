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
    bool m_bits[10000000];
    uint64_t b_size;

    public:
        BloomFilter(uint64_t capacity,double errorRate);
        void initBloom(uint64_t capacity,double errorRate);
        void add(string data);
        void printSetBits();
        void    setBSize(uint64_t size);
        bool possiblyContains(string data) const;
        uint64_t getMBitsSize();
        uint64_t getNumHashes();
        uint64_t getBSize();
        uint64_t hash(string data,uint64_t filterSize,int n) const;
        ~BloomFilter();

};

#endif
