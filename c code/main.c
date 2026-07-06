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

#define FILE "learn_mouse"
#define FRECS "_freqs"
#define ENTROPIS "_entropies"
#define N_LINES 373

int main() {

    /*Initializations*/
    initialize_properties_matrix();

    Chain ch = file_megaAacids(SEQS FILE TXT , N_LINES);

    print_chain_to_file(ch, RESULTS FILE FRECS TXT);
    
    Entropies *entropy = malloc(CHAINLEN * sizeof(Entropies));
    all_entropies(ch, entropy, .5);

    print_entropies_to_file(entropy, RESULTS FILE ENTROPIS TXT);

    /*Free memory*/
    free_properties_matrix();
    return 0;
}