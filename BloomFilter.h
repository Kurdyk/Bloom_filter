#ifndef BLOOM_FILTERS_BLOOMFILTER_H
#define BLOOM_FILTERS_BLOOMFILTER_H

#include <iostream>
#include <vector>

using namespace std;


class BloomFilter {

    private:uint64_t size;
    private:uint nf;
    private:std::vector<bool>* tab;


    public:BloomFilter(uint64_t n, uint nf) {
        this->size = n;
        this->nf = nf;
        tab = new vector<bool>(n, false);
    }

    public:void add_element(uint64_t value) {
        uint64_t hashes[this->nf];
        memset(hashes, 0, nf);
        multihash(value,  hashes);
        for (int i = 0; i < this->nf; i++) {
            this->tab->at(hashes[i]) = true;
        }
    }

    public:bool is_present(uint64_t value) {
        uint64_t hashes[this->nf];
        multihash(value,  hashes);
        for (int i = 0; i < this->nf; i++) {
            if(!this->tab->at(hashes[i])) return false;
        }
        return true;
    }

    /** Hash function that uses xor properties
    * @param x a 64-bits value to hash
    * @return hashed value
    */
    private:uint64_t xorshift64(uint64_t x) {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        return x;
    }

    /** Generate multiple hashes of the same value using sequences of xorshifts.
    * @param x The value to hash
    * @param hashes An array already allocated to fill with hash values
    */
    private:void multihash(uint64_t x, uint64_t * hashes) {
        // TODO : gérer la valeur 0

        // Init 64 bits hashes
        hashes[0] = xorshift64(x);
        for (uint64_t i=1 ; i < this->nf ; i++)
            hashes[i] = xorshift64(hashes[i-1]);

        for (uint64_t i=0 ; i < this->nf ; i++)
            hashes[i] %= this->size;
    }


    public:void print() {
        for (uint64_t i = 0; i < this->size; i++) {
            cout << tab->at(i) << ' ';
        }
    }

    public:void delete_tab() {
        delete this->tab;
    }



};


#endif //BLOOM_FILTERS_BLOOMFILTER_H
