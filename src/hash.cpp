#include "BloomFilter.h"
#include "fastaProcessing.h"


using namespace std;

int main(int argc, char *argv[]) {

    if (argc != 6) {
        cout << "Usage : ./hash file_name k n nf r \n"
                "Where :\n -file_name is the path to the fasta file.\n"
                " -k is the length of the k-mers.\n"
                " -n is the size in bits of the BloomFilter.\n"
                " -nf is the number of hash done.\n"
                " -r is the number of random request to do." << endl;
        return 0;
    }

    /// Getting command line arguments and setting the global variables.
    int k = atol(argv[2]);
    uint64_t n = atol(argv[3]);
    int nf = atoi(argv[4]);
    uint64_t r = atol(argv[5]);


    FILE* file = fopen(argv[1], "r");
    BloomFilter BF = BloomFilter(n, nf);

    process_fasta(file, BF, k);
    fclose(file);

    random_requests(r, BF);


    return 0;
}
