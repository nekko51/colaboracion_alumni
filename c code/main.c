#include "head.h"

/*  lines in each seq file:
learn human:    1309
learn mouse:    373
test human:     1388
test mouse:     1379
*/

#define SEQS "seqs/"
#define RESULTS "results/"
#define TXT ".txt"

#define FILE "learn_human"
#define MEGAAA "_freqs"
#define N_LINES 1309

int main() {

    Chain ch = file_megaAacids(SEQS FILE TXT , N_LINES);

    print_chain_to_file(ch, RESULTS FILE MEGAAA TXT);
    
    Entropies *entropy = malloc(CHAINLEN * sizeof(Entropies));
    all_entropies(ch, entropy, .5);

    return 0;
}