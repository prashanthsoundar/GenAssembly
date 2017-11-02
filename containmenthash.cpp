#include<iostream>
#include<climits>
#include<vector>
#include <string>
#include "bloomfilter.hpp"
#include "ntHashIterator.hpp"
using namespace std;
class ContainmentHash{
	void printhashes(const size_t hashes[]){
		for(size_t i = 0; i <100;i++){
			cout<<"---"<<hashes[i]%1000<<endl;
		}
	}
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
	
		ntHashIterator kmerhashes(se1,number_of_hashes,kmer_size);	
		int count = 0;
		while (kmerhashes != kmerhashes.end()) {
			cout<<bloom.contains(*itr1)<<endl;
			printhashes(*kmerhashes);
						
			if(bloom.contains(*kmerhashes)){
				count++;
			}
			++itr1;
		}
		
		cout<<"Count:"<<count;
	}
	
};
