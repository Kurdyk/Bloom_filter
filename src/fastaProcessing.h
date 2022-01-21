#ifndef BLOOM_FILTERS_FASTAPROCESSING_H
#define BLOOM_FILTERS_FASTAPROCESSING_H

#include <iostream>
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <bitset>

#include "BloomFilter.h"

void process_fasta(FILE* fasta, BloomFilter &BF, const int &length);
void random_requests(const uint64_t &nb_requests, BloomFilter &BF);


#endif //BLOOM_FILTERS_FASTAPROCESSING_H
