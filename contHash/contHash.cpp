#include<iostream>
#include "include/BloomFilter.h"
#include "include/MinhashNaive.h"
#include<math.h>
#include<string>
typedef unsigned long long int uint64_t;
using namespace std;

//vars for minHash count estimator
uint64_t prime = 9999999999971UL;
uint64_t kSize = 7;
//uint64_t h = 10;
double p = 0.01;


//to be initiliased
string l = "CATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAGCATGGACCGACCAG";
string s = "GCAGTACCGATCGT";


int main()
{
    //cin>>largeString;
    BloomFilter *b = new BloomFilter(1.15*l.length(),p);
    MinhashNaive *m = new MinhashNaive(b->getNumHashes(),b->getMBitsSize());
    
    vector<string> kmers = m->computeHashedKmers(s,kSize);
    
    uint64_t bSize=0;
    for(int i=0;i<l.length()-kSize+1;i++)
    {
        //this loop is slow.. Make it faster!
        if(!b->possiblyContains(l.substr(i,kSize)))
        {
            b->add(l.substr(i,kSize));
            bSize++;
        }
        
    }
   cout<<bSize<<endl;
    uint64_t intersectionCount=0;
    for(int i=0;i<kmers.size();i++)
    {
        if(b->possiblyContains(kmers[i]))
        {
            intersectionCount++;
        }
    }
    
    intersectionCount-=(uint64_t)floor(p*b->getNumHashes());
    cout<<intersectionCount<<endl;
    
    
    
    
    
    
    
    //construct a bloom filter //
    
    
   

    return 0;
}
