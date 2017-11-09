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
uint64_t kSize = 20;
uint64_t h = 200;
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
float jaccardindex( vector<uint64_t> &v1, vector<uint64_t> &v2)
{

	vector<uint64_t> v3;
	vector<uint64_t> v4;
	int sizeintersection;
	int sizeunion;

	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());

	set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
	std::set_union (v1.begin(),v1.end(),v2.begin(),v2.end(), back_inserter(v4));
	sizeintersection = v3.size();
	sizeunion = v4.size();

	float ji = ((float)sizeintersection) / sizeunion;

	return ji;

}
/*int main(){
  int number_of_hashes[] = {1,5,10,100,1000,1100,1200,1500,2000,2500};
  int kmer_size = 20;
  for(int i:number_of_hashes){
  h = i;
  containHash();	
  }
  }*/
int containHash()
{
	/*	ifstream input("input.txt");
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
	 */
	BloomFilter *b = new BloomFilter(1.15*l.length(),p);
	MinhashNaive *m1 = new MinhashNaive(h,9999999999971UL);
	vector<string> kmers = m1->computeHashedKmers(s,kSize); 
	//vector<uint64_t> sk1 = m1->getSketch();
	// MinhashNaive *m2 = new MinhashNaive(h,9999999999971UL);
	// vector<string> kmers2 = m2->computeHashedKmers(s,kSize);
	//vector<uint64_t> sk2 = m2->getSketch();
	//    cout<<"======"<<sk1<<endl<<sk2<<endl;
	//for(int y=0;y<sk1.size();y++)
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

		cout<<"String1 Size:"<<l.length()<<" "<<bSize<<endl;
		cout<<"String2 Size:"<<s.length()<<" "<<sizeSmallStr<<endl;
	//	cout<<"K-Mer Size\t\t: "<<kSize<<endl;
	//	cout<<"Hash Count\t\t: "<<h<<endl;
	intersectionCount-=(uint64_t)floor(p*h);
	//	cout<<"Intersection Est.\t: "<<intersectionCount<<endl;

	float containmentEst=intersectionCount/float(h);
	//	cout<<"Containment Est.\t: "<<containmentEst<<endl;

	float jaccardEst = (sizeSmallStr * containmentEst) / (sizeSmallStr + bSize - sizeSmallStr * containmentEst);
	cout<<"Jaccard Est. by containmentHash: "<<jaccardEst<<endl;
	float trueJaccard = jaccardIndex(l,s);
	cout<<"Actual Jaccard Est.: "<<trueJaccard<<endl;
	float relError = abs(jaccardEst-trueJaccard)/trueJaccard;
	cout<<"Relative Error of containment hash: "<<setprecision(6)<<relError<<endl;
	cout<<"------------------------------------"<<endl;
	//    ofstream in("output.csv");

	return 0;
}
int main(int argc,char* argv[]){
	if(argc>1){
		kSize = atoi(argv[1]);
		h = atoi(argv[2]);
	}
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
	MinhashNaive *m1 = new MinhashNaive(h,9999999999971UL);
	m1->computeHashedKmers(l,kSize);
	vector<uint64_t> sk1 = m1->getSketch();
	MinhashNaive *m2 = new MinhashNaive(h,9999999999971UL);
	m2->computeHashedKmers(s,kSize); 
	vector<uint64_t> sk2 = m2->getSketch();
	cout<<"KmerSize: "<<kSize<<endl<<"Number of Hash Values: "<<h<<endl<<"Minhash similarity: "<<::jaccardindex(sk1,sk2)<<endl;
	containHash();

}
