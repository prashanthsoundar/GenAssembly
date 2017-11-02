
#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <exception>
#include <string>
#include <iostream>
#include <climits>
using namespace std;

#include "MurmurHash3.h"
#include "containmenthash.cpp"

uint64_t generate_hash(const string& kmer, const uint32_t seed) {
    uint64_t out[2];
    MurmurHash3_x64_128((void *)kmer.c_str(), kmer.size(), seed, &out);
    return out[0];
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


class KminHash {
public:
    int size;
    int seed;
    int num_hash;
    std::map<std::string, vector<uint64_t> > m;
    std::vector<uint64_t> sketch;

    KminHash(int k, int h) {
        size = k;
        seed = 0;
        num_hash = h;

        for (int i = 0; i < num_hash; i++) {
            sketch.push_back(ULONG_MAX);
        }
    }

    // TODO: Check if its working fine
    static uint64_t xorShift64(uint64_t a) {
        a ^= (a << 21);
        a ^= (a >> 35);
        a ^= (a << 4);
        return a;
    }
    
    //Ref : https://en.wikipedia.org/wiki/Xorshift
    uint64_t xorshift64star(uint64_t x) {
        //uint64_t x = state[0];
        x ^= x >> 12; // a
        x ^= x << 25; // b
        x ^= x >> 27; // c
        //state[0] = x;
        return x * 0x2545F4914F6CDD1D;
    }

    void add_kmer(const std::string& kmer) {
        vector<uint64_t> kmer_array;
        uint64_t temp = generate_hash(kmer, seed);
        kmer_array.push_back(temp);
        if (sketch[0] > temp) {
            sketch[0] = temp;
        }

        for (int i = 1; i < num_hash; i++) {
            temp = xorshift64star(kmer_array[i-1]);
            kmer_array.push_back(temp);
            if (sketch[i] > temp) {
                sketch[i] = temp;
            }
        }
        m[kmer] = kmer_array;
    }

    void generate_kmer(const string s) {
        for (int i = 0; i <= s.length() - size; i++) {
            const string sub_string = s.substr(i, size);
            add_kmer(sub_string);
        }
    }

    void compare(KminHash* other) {
        cout << "Sketch1" << endl;
        for (int i = 0; i < num_hash; i++) {
            cout << sketch[i] << ",";
        }
        cout << endl;
        cout << "Sketch2" << endl;
        for (int i = 0; i < num_hash; i++) {
            cout << other->sketch[i] << ",";
        }
        cout << endl;
        cout << ::jaccardindex(sketch, other->sketch) << endl;
    }
};


int main() {

    KminHash* temp1 = new KminHash(3, 17);
    temp1->generate_kmer("CATGGACCGACCAG");

    KminHash* temp2 = new KminHash(3, 17);
    temp2->generate_kmer("GCAGTACCGATCGT");

    temp1->compare(temp2);
    ContainmentHash* cmh = new ContainmentHash();
    cmh->setvalues(3,100);
    cmh->calculatesimilarity("CATGGACCGACCAG","GCAGTACCGATCGT");
}
