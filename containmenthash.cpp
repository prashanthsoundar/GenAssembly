#include<iostream>
#include<climits>
#include<vector>
#include <string>
#include "bloomfilter.hpp"
#include "ntHashIterator.hpp"
using namespace std;
class ContainmentHash{
	
	int kmer_size;
	int number_of_hashes;
public:
	void setvalues(int size,int no_hash){
		kmer_size = size;
		number_of_hashes = no_hash;
	}
	float calculatesimilarity(const string genome1,const string genome2){
		BloomFilter bloom(1000,number_of_hashes,kmer_size);	
		string se,se1;
		if(genome1.length()>genome2.length()){
			se = genome1;
			se1 = genome2;
		}else{
			se = genome2;
			se1 = genome1;
		}
		ntHashIterator itr(se,number_of_hashes, kmer_size);
		while (itr != itr.end()) {
			bloom.insert(*itr);
			++itr;
		}
	
//		bloom.insert("ddd");
//		bloom.storeFilter("filter.bf");	
//		BloomFilter bloom("filter.bf");
		ntHashIterator itr1(se1,number_of_hashes, kmer_size);
		int count = 0;
		while (itr1 != itr1.end()) {
			cout<<bloom.contains(*itr1)<<endl;
			if(bloom.contains(*itr1)){
				count++;
			}
			++itr1;
		}
		
		cout<<"Count:"<<count;
	}
	
};
