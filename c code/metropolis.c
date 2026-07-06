#include "head.h"

double hamming_distance(const char* seq1, const char* seq2, int n) {
    double dist = 0.0;
    for(int i=0; i<n; i++) {
        if (seq1[i] != seq2[i]) dist++;
    }
    return(dist);
}

double calculate_energy(const char* seq1, const char* seq2, int n) {
    return(hamming_distance(seq1, seq2, n));
}

/*We're probably interested in making selective sweeps, not sweeps of the whole sequence, but this'll do for now*/
/*Do we want a sweep of the whole sequence, or random position sweeps?*/
void metropolis_sweep(char* murine_seq, const char* human_ref, double beta, int n) {
    for(int i=0; i<n; i++) {
        double old_energy = calculate_energy(murine_seq, human_ref, n);

        /*We propose a tentative change:*/
        char old_aa = murine_seq[i];
        char new_aa;
        do {
            int idx = 1*fran*1;
            new_aa = AMINOACIDS[idx]; //Lógicamente idx ha de ser un entero; depende de qué rand usemos (parisi-rapuano, etc, habrá q cambiarlo)
        } while(new_aa == old_aa);

        /*Metropolis acceptance*/
        double new_energy = calculate_energy(murine_seq, human_ref, n);//As with Ising, this could be way more efficient since there's a finite number of delta_e's
        double delta_e = new_energy - old_energy;
        if (delta_e < 0 || fran < exp(-beta*delta_e)) {
            murine_seq[i] = new_aa;
            old_energy = new_energy;
        }
    }

}