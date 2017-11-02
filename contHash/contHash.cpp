#include<iostream>
#include "include/BloomFilter.h"
#include<string>
typedef unsigned long long int uint64_t;
using namespace std;

//vars for minHash count estimator
uint64_t prime = 9999999999971UL;
uint64_t kSize = 11;
uint64_t h = 10;
double p = 0.01;


//to be initiliased
string smallString;
string largeString;


int main()
{
    //cin>>largeString;
    BloomFilter *b = new BloomFilter(300000,p);
    b->add("testWord");
    b->printSetBits();
    //construct a bloom filter //
    
    
   

    return 0;
}
