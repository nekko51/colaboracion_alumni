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

    /*Variable Initializations*/
    initialize_properties_matrix();

    /*Code preparation*/
    Chain ch = file_megaAacids(SEQS FILE TXT , N_LINES);
    if(ch.aas->elements[0] < -0.5) {
        fprintf(stderr, "Error: Chain returned negative values (%lf); stopping execution...\n", ch.aas->elements[0]);
        return 1;
    }

    print_chain_to_file(ch, RESULTS FILE FRECS TXT);
    
    Entropies *entropy = malloc(CHAINLEN * sizeof(Entropies));
    all_entropies(ch, entropy, .5);

    print_entropies_to_file(entropy, RESULTS FILE ENTROPIS TXT);

    /*Free memory*/
    
    return 0;
}