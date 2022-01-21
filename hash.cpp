#include <iostream>
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <bitset>

#include "BloomFilter.h"

using namespace std;

/** Global variables declaration.
 * encoding is a map to give value to each nucleotide.
 * rev_table is a map that associate a nucleotide with it's complement
 * k is the length of the wanted k-mers
 * filter is used to make bits to bits & necessary to calculate hash values from precedent value.
 */
unordered_map<char, uint64_t> encoding;
unordered_map<char, char> rev_table;
int k;
uint64_t filter;

/**
 * Set the file pointer to the start of the DNA sequence.
 * @param fasta a file in fasta format.
 */
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

/**
 * Read the next character in a fasta file excluding end of line and 'N' characters.
 * @param fasta a file in fasta format.
 * @return the next character of the sequence.
 */
char read_fasta(FILE* fasta) {
    int c = fgetc(fasta);
    if (c == '\n' || c == 'N') return read_fasta(fasta);
    return (char) c;
}

/**
 * Initialise the filter.
 * @return the value of the filter linked to the size k.
 */
uint64_t set_filter() {
    string str;
    str.resize(2*(k+1), '1');
    std::fill(str.begin(), str.begin() + 2, '0');
    bitset<64> filtre(str);
    return filtre.to_ullong();
}

/**
 * Choose between the to form of a k-mer between itself and it's reverse complement.
 * @param kmer the actual kmer.
 * @return the smallest in lexicographic order between kmer and it's reverse complement.
 */
string choose_complement(const string &kmer) {
    string reverse;
    for (int i = 0; i < k; i++) {
        reverse += rev_table.at(kmer.at(k - i - 1));
    }
    return (reverse < kmer) ? reverse : kmer;
}

/**
 * Hash the given kmer.
 * @param kmer a kmer to hash.
 * @return the minimum between hash of kmer or of it's reverse complement.
 */
uint64_t hash_string(string kmer) {
    uint64_t value = 0;
    kmer = choose_complement(kmer);
    for (int i = 0; i < k; i++) {
        value = (value << 2) + encoding.at(kmer[i]);
    }
    return value;
}

/**
 * Hash the reverse complement of a kmer.
 * @param kmer a kmer to reverse and hash.
 * @return the hash of kmer's complement.
 */
uint64_t hash_rev(string kmer) {
    string reverse;
    for (int i = 0; i < k; i++) {
        reverse += rev_table.at(kmer.at(k - i - 1));
    }
    uint64_t value = 0;
    for (int i = 0; i < k; i++) {
        value = (value << 2) + encoding.at(reverse[i]);
    }
    return value;
}

/**
 * Hash the first k-mer of the fasta file and it's complement for initialisation.
 * @param kmer the first k-mer of the file.
 * @param prev_val a reference to the hash of the previous k-mer (without considering it's complement).
 * @param prev_val_rc a reference to the hash of the previous k-mer's complement.
 * @return the hash of the first k-mer of the file (or it's complement).
 */
uint64_t hash_start(string kmer, uint64_t &prev_val, uint64_t &prev_val_rc) {
    uint64_t value = 0;
    for (int i = 0; i < k; i++) {
        value = (value << 2) + encoding.at(kmer[i]);
    }
    prev_val = value;
    prev_val_rc = hash_rev(kmer);
    return hash_string(kmer);
}

/**
 * Determine the hash of the actual k-mer from the previous k-mer, it's hash value,
 * the hash value of it's complement and the next nucleotide of the sequence.
 * @param prev_val a reference to the hash of the previous k-mer (without considering it's complement).
 * @param prev_val_rc a reference to the hash of the previous k-mer's complement.
 * @param prev_kmer a reference to the previous k-mer in the file.
 * @param new_end the nucleotide just after prev_kmer determining the actual k-mer.
 * @return the hash value of the actual k-mer.
 */
uint64_t hash_from_previous_kmer(uint64_t &prev_val, uint64_t &prev_val_rc, string &prev_kmer, const char &new_end){
    prev_kmer.erase(0, 1).push_back(new_end);
    prev_val = ((prev_val << 2) & filter) + encoding.at(new_end);
    prev_val_rc = (((prev_val_rc >> 2) + (encoding.at((rev_table.at(new_end))) << (2 * (k - 1)))));
    return (prev_val < prev_val_rc)? prev_val : prev_val_rc;
    /*
    if (choose_complement(prev_kmer) == prev_kmer) {
        return prev_val;
    }
    else {
        return prev_val_rc;
    }
     */
}

/**
 * Hash and add to a BloomFilter every k-mer of a fasta file.
 * @param fasta a file in the fasta format.
 * @param BF a BloomFilter to store the k-mer.
 */
void process_fasta(FILE* fasta, BloomFilter &BF) {
    cout << "Running..." << endl;
    filter = set_filter();
    set_fasta(fasta);

    long n = 0;
    char c;
    string kmer;
    bool wanted_size = false;
    uint64_t value;
    uint64_t prev_value;
    uint64_t prev_value_rc;
    while((c = read_fasta(fasta)) != EOF) {
        if(!wanted_size) {
            kmer.push_back(c);
            if (kmer.length() == k) {
                wanted_size = true;
                value = hash_start(kmer, prev_value, prev_value_rc);
                BF.add_element(value);
                n++;
            }
        } else {
            value = hash_from_previous_kmer(prev_value, prev_value_rc, kmer, c);
            BF.add_element(value);
            n++;
        }
    }
    cout << "End of file. " << n << " elements added."<< endl;
}

/**
 * Generate a random k-mer.
 * @return a string containing a random k-mer.
 */
string random_kmer() {
    string kmer;
    char possible_value[] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < k; i ++) {
        kmer += possible_value[rand() % 4];
    }
    return kmer;
}

/**
 * Make nb_requests of random k-mer to a BloomFilter.
 * @param nb_requests the wanted number of requests.
 * @param BF a BloomFilter.
 */
void random_requests(const uint64_t &nb_requests, BloomFilter &BF){
    for (uint64_t i = 0; i < nb_requests; i++) {
        string to_test = random_kmer();
        cout << "test if kmer (or it's reverse complement) : "
        << to_test << " is present : " << BF.is_present(hash_string(to_test)) << endl;
    }
}


int main(int argc, char *argv[]) {

    /// Getting command line arguments and setting the global variables.
    k = atol(argv[2]);
    uint64_t n = atol(argv[3]);
    int nf = atoi(argv[4]);
    uint64_t r = atol(argv[5]);


    encoding['A'] = 0; // 00
    encoding['T'] = 3; // 11
    encoding['C'] = 1; // 01
    encoding['G'] = 2; // 10

    rev_table['A'] = 'T';
    rev_table['T'] = 'A';
    rev_table['C'] = 'G';
    rev_table['G'] = 'C';


    FILE* file = fopen(argv[1], "r");
    BloomFilter BF = BloomFilter(n, nf);

    process_fasta(file, BF);
    fclose(file);

    srand(time(NULL));
    random_requests(r, BF);


    return 0;
}
