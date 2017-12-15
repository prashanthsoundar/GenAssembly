//
// Created by Naga Srinikhil Reddy on 12/14/17.
//
#include<iostream>
#include "include/BloomFilter.h"
#include "include/MinhashNaive.h"
#include<math.h>
#include<set>
#include <fstream>
#include<string>
#include<iomanip>
#include<set>
#include<stdlib.h>
#include <iostream>
#include <fstream>

BloomFilter *b;
int h = 200;
int kSize = 20;
double p = 0.001;
vector <string> r1;
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
void findSim(string s){

    MinhashNaive *m1 = new MinhashNaive(h,9999999999971UL);
    vector<string> kmers = m1->computeHashedKmers(s,kSize);
    uint64_t intersectionCount = 0;
    set<string> temp;
    for(int i =0;i<s.length()-kSize+1;i++){
        temp.insert(s.substr(i,kSize));
    }
    uint64_t sizeSmallStr=temp.size();
    for(int i=0;i<kmers.size();i++)
    {
        if(b->possiblyContains(kmers[i]))
        {
            intersectionCount++;
        }
    }
    intersectionCount-=(uint64_t)floor(p*h);
    float containmentEst=intersectionCount/float(h);
    uint64_t bSize = b->getBSize();
    float jaccardEst = (sizeSmallStr * containmentEst) / (sizeSmallStr + bSize - sizeSmallStr * containmentEst);
    cout<<"Jaccard Est. by containmentHash: "<<jaccardEst<<endl;
    cout<<"------------------------------------"<<endl;
}
int main(int argc,char* argv[]){
    string file_input;

    if(argc>1) {
        file_input = argv[1];
    }
    ifstream input(file_input);
    if(!input.is_open()){
        cout<<"File name not correct. Exiting"<<endl;
        return 0;
    }
    string f_l,line;
    int q = 0;
    while( std::getline( input, line ) ){
        if(q == 0){
            q=1;
            continue;
        }
        if( line.empty() || line[0] == '>' ){
            r1.push_back(f_l);
            f_l = "";
        }
        else{
            f_l+=line;
        }
    }
    cout<<r1.size();
    ifstream file_obj;
    file_obj.open("bloom.obj", ios::in);

    b =  (BloomFilter*)malloc(sizeof(BloomFilter));
    file_obj.read((char*)b, sizeof(BloomFilter));
//    cout<<"Number of bits"<<b->getBSize()<<endl;
  //  cout<<"Number of bits"<<b->getMBitsSize()<<endl;
  //  cout<<"Number of hashes"<<b->getNumHashes()<<endl;
    float  t;
   cout<<"Bloom Filter loaded from file"<<endl; 
   for (int k = 0; k < r1.size(); ++k) {
        std::time_t t1 = std::time(nullptr);
        findSim(r1.at(k));
        std::time_t t2 = std::time(nullptr);
        t += (t2-t1);
    }
//    cout<<"Total run time for "<<r1.size()<<" reads each of size "<<r1.at(0).size()<<" is "<<(t)<<endl;

}
