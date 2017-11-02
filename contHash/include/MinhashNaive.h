#ifndef _minhashnaive_h
#define _minhashnaive_h
#include<vector>
#include<string>

using namespace std;

typedef unsigned long long int uint64_t;


class MinhashNaive
{
    uint64_t numHashes;
    vector<uint64_t> sketch;
    vector<string> sketchStr;
    vector<string> hashedKmers;
    uint64_t largeMod;
    
    public:
        MinhashNaive(uint64_t numHashes,uint64_t largeMod);
        vector<string> computeHashedKmers(string s,uint64_t kmerSize);
        uint64_t hash(string data,int n) const;
        ~MinhashNaive();
    
    
};


#endif
