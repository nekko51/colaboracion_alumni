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
void initialize(Chain* human_ref) {
    ini_ran(time(NULL));
    initialize_properties_matrix();
    initialize_char_to_int_LUT();
    file_megaAacids(SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES, human_ref);
}

int main() {//i'm so happy i don't have to free every single malloc'd array if there's an error; if we had to, I think it's the only use of GOTO that wouldn't get you fired
    /*Initial parameters*/
    int n_betas = 50;
    int n_sweeps = 800;
    int n_metropolis = 20;
    int n_entropies = CHAINLEN;
    double entropy_order_q = 0.34;
    double scale_factor = 1.0;
    double cooling_rate = 1.02;
    double weighs[8] = {0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0};
    Chain human_ref;

    /*Variable & RNG Initializations*/
    initialize(&human_ref);
    
    double** betas = malloc(n_betas*sizeof(double*));
    if(betas == NULL) {
        fprintf(stderr, "Couldn't assign memory to betas matrix\n");
        return(1);
    }
    for(int i=0; i<n_betas; i++) {
        betas[i] = malloc(CHAINLEN*sizeof(double));
        if(betas[i] == NULL) {
            fprintf(stderr, "Couldn't assign memory to betas column %d/%d\n", i+1, n_betas);
            return(1);
        }
    }

    Entropies* oll_entropies = malloc(CHAINLEN*sizeof(Entropies));
    if(oll_entropies == NULL) {
        fprintf(stderr, "Couldn't assign memory to all_entropies array\n");
        return(1);
    }
    double* entropies = malloc(CHAINLEN*sizeof(double));
    if(entropies == NULL) {
        fprintf(stderr, "Couldn't assign memory to entropies array\n");
        return(1);
    }

    all_entropies(&human_ref, oll_entropies, entropy_order_q);
    weigh_entropies(oll_entropies, entropies, weighs);//unfinished function, currently returns 1/2*(saa+spp)
    print_entropies(oll_entropies);
    generate_betas(betas, n_betas, entropies, n_entropies, scale_factor, EPSILON, cooling_rate);
    mega_metropolis(SEQS FILE_L_MOUSE TXT, SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES, n_sweeps, betas, n_betas, n_metropolis);




    /*Code preparation*/
    // Chain ch = file_megaAacids(SEQS FILE_L_MOUSE TXT, L_MOUSE_N_LINES);
    // if(negative_chain(&ch) == 1) return 1;

    // print_chain_to_file(ch, RESULTS FILE_L_MOUSE FRECS TXT);
    
    // Entropies *entropy = malloc(CHAINLEN * sizeof(Entropies));
    // all_entropies(ch, entropy, .5);

    // print_entropies_to_file(entropy, RESULTS FILE_L_MOUSE ENTROPIS TXT);

    /*Free memory*/
    for(int i=0; i<n_betas; i++) {
        free(betas[i]);
    }
    free(betas);
    free(oll_entropies);
    free(entropies);
    
    return 0;
}