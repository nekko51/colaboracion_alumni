#include "head.h"

/*  lines in each seq file:
learn human:    1309
learn mouse:    373
test human:     1388
test mouse:     1379
*/

/** 
 * @brief file definitions in head.h
 * 
 * @param SEQS "seqs/"
 * @param RESULTS "results/"
 * @param TXT
 * @param FILE_L_MOUSE "learn_mouse"
 * @param L_MOUSE_N_LINES 373
 * @param FILE_L_HUMAN "learn_human"
 * @param L_HUMAN_N_LINES 1309
 * @param FRECS "_freqs"
 * @param ENTROPIS "_entropies"
 *
*/

//returns 1 if there's an error, 0 if everything went smoothly
int run_metropolis() {
    Chain human_ref_seq = file_megaAacids(SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES);
    if(negative_chain(human_ref_seq) == 1) return(1);

    Chain murine_seed = generate_murine_seed(megachain);
    if(negative_chain(murine_seed) == 1) return(1);
}

int main() {

    /*Variable Initializations*/
    initialize_properties_matrix();

    /*Code preparation*/
    Chain ch = file_megaAacids(SEQS FILE_L_MOUSE TXT, L_MOUSE_N_LINES);
    if(negative_chain(ch) == 1) return 1;

    print_chain_to_file(ch, RESULTS FILE_L_MOUSE FRECS TXT);
    
    Entropies *entropy = malloc(CHAINLEN * sizeof(Entropies));
    all_entropies(ch, entropy, .5);

    print_entropies_to_file(entropy, RESULTS FILE_L_MOUSE ENTROPIS TXT);

    /*Free memory*/
    
    return 0;
}