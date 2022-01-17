#include <iostream>
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <cstring>
#include <bitset>

#include "BloomFilter.h"

using namespace std;

void set_fasta(FILE* fasta){
    int c;
    c = fgetc(fasta);
    while(c != EOF && c != '\n') {
        c = fgetc(fasta);
    }
    if (c == EOF) {
        perror("Bad file format");
        exit(1);
    }
}


char read_fasta(FILE* fasta) {
    int c = fgetc(fasta);
    if (c == '\n' || c == 'N') return read_fasta(fasta);
    return (char) c;
}


uint64_t set_filter(int k) {
    std::string str;
    str.resize(2*(k+1), '1');
    std::fill(str.begin(), str.begin() + 2, '0');
    std::bitset<64> filter(str);
    return filter.to_ullong();
}


uint64_t quick_pow(uint64_t base, uint64_t power)
{
    uint64_t value = 1;
    while (power)
    {
        if (power % 2)
            value *= base;
        power /= 2;
        base *= base;
    }
    return value;
}

uint64_t hash_string(const string str, int length , unordered_map<char, int> encoding) {
    uint64_t value = 0;
    for (int i = 0; i < length; i++) {
        value += quick_pow(4, length - 1 - i) * encoding.at(str[i]);
    }
    return value;
}


u_int64_t hash_from_previous_kmer(uint64_t prev_val, char new_end, uint64_t filter, unordered_map<char, int> encoding){
    return ((prev_val << 2) & filter) + encoding.at(new_end);
}


void process_fasta(FILE* fasta, int k, BloomFilter BF,  const unordered_map<char, int>& encoding) {

    uint64_t filter = set_filter(k);

    set_fasta(fasta);

    char c;
    std::string kmer;
    bool wanted_size = false;
    uint64_t value;
    while((c = read_fasta(fasta)) != EOF) {
        if(!wanted_size) {
            kmer.push_back(c);
            if (kmer.length() == k) {
                wanted_size = true;
                // TODO : choose the good kmer.
                value = hash_string(kmer, k, encoding);
                BF.add_element(value);
            }
        } else {
            // TODO : choose the good kmer and shit with the new c.
            value = hash_from_previous_kmer(value, c, filter, encoding);
            BF.add_element(value);
        }
    }
}


std::string random_kmer(int k) {
    std::string kmer;
    char possible_value[] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < k; i ++) {
        kmer += possible_value[rand() % 4];
    }
    return kmer;
}


void random_requests(uint64_t nb_requests, BloomFilter BF, int k,  const unordered_map<char, int>& encoding){
    for (uint64_t i = 0; i < nb_requests; i++) {
        cout << BF.is_present(hash_string(random_kmer(k), k, encoding)) << endl;
    }
}

int main(int argc, char *argv[]) {
    int k = 10;
    srand(time(NULL));

    unordered_map<char, int> encoding;
    encoding['A'] = 0;
    encoding['C'] = 1;
    encoding['T'] = 2;
    encoding['G'] = 3;


    FILE* file = fopen(argv[1], "r");
    BloomFilter BF = BloomFilter(250000, 12);

    process_fasta(file, k, BF, encoding);

    fclose(file);
    BF.delete_tab();

    return 0;
}
