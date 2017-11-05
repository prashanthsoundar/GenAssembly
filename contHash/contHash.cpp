#include<iostream>
#include "include/BloomFilter.h"
#include "include/MinhashNaive.h"
#include<math.h>
#include<string>
#include<iomanip>
#include<set>
#include<stdlib.h>
#include <iostream>
#include <fstream>
typedef unsigned long long int uint64_t;
using namespace std;

//vars for minHash count estimator
uint64_t prime = 9999999999971UL;
uint64_t kSize = 15;
uint64_t h = 100;
double p = 0.001;
string l,s,line;
float jaccardIndex(string a,string b)
{
    set<string> s1;
    set<string> s2;
    std::pair<std::set<string>::iterator,bool> ret;
    uint64_t unionCount=0;
    uint64_t dupCount=0;
    
    //remove duplicates from a
    for(int i=0;i<a.length()-kSize+1;i++)
    {
        s1.insert(a.substr(i,kSize));
    }
    
    s2=s1;
    s1.clear();
    
    // remove duplicates from b
    for(int i=0;i<b.length()-kSize+1;i++)
    {
        ret = s1.insert(b.substr(i,kSize));
        if(!ret.second)
            dupCount++; //this may contain duplicates of b string
    }
    
    s1.clear(); //this is not absolutely necessary
    
    for(int i=0;i<b.length()-kSize+1;i++)
    {
        ret = s2.insert(b.substr(i,kSize));
        if(!ret.second)
            unionCount++;
    }
    
    unionCount-=dupCount; //actual UnionCount
    
    float ji = (float)unionCount/s2.size();
    return ji;
    
}

int main()
{
    	ifstream input("input.txt");
        while( std::getline( input, line ) ){
        	if( line.empty() || line[0] == '>' ){
			 continue;
                }   
                else{
                        if(l.empty()){
                                l = line;
                        }   
                        else{
                                s = line;
                        }   
                }   
        }   
    BloomFilter *b = new BloomFilter(1.15*l.length(),p);
    MinhashNaive *m = new MinhashNaive(h,9999999999971UL);
    vector<string> kmers = m->computeHashedKmers(s,kSize); 
    set<string> temp;
    for(int i =0;i<s.length()-kSize+1;i++){
        temp.insert(s.substr(i,kSize));
    }
    uint64_t sizeSmallStr=temp.size();
    temp.clear();
    for(int i=0;i<l.length()-kSize+1;i++){
            b->add(l.substr(i,kSize));
            temp.insert(l.substr(i,kSize));
    }
    uint64_t bSize=temp.size();
    temp.clear();
    
    uint64_t intersectionCount=0;
    for(int i=0;i<kmers.size();i++)
    {
        if(b->possiblyContains(kmers[i]))
        {
            intersectionCount++;
        }
    }
    
    cout<<"K-Mer Size\t\t: "<<kSize<<endl;
    cout<<"Hash Count\t\t: "<<h<<endl;
    intersectionCount-=(uint64_t)floor(p*h);
    cout<<"Intersection Est.\t: "<<intersectionCount<<endl;
   
    float containmentEst=intersectionCount/float(h);
    cout<<"Containment Est.\t: "<<containmentEst<<endl;
    
    float jaccardEst = (sizeSmallStr * containmentEst) / (sizeSmallStr + bSize - sizeSmallStr * containmentEst);
    cout<<"Jaccard Est. by containmentHash: "<<jaccardEst<<endl;
    float trueJaccard = jaccardIndex(l,s);
    cout<<"Actual Jaccard Est.\t: "<<trueJaccard<<endl;
    float relError = abs(jaccardEst-trueJaccard)/trueJaccard;
    cout<<"Relative Error\t\t: "<<setprecision(6)<<relError<<endl;
    cout<<"------------------------------------"<<endl;
    return 0;
}
