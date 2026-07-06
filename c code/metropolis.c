#include "head.h"

#define fran 3

double hamming_distance(const char* seq1, const char* seq2, int n) {
    double dist = 0.0;
    for(int i=0; i<n; i++) {
        if (seq1[i] != seq2[i]) dist++;
    }
    return(dist);
}

//Calculates the humanness of a sequence based on the schrödinger human AA chain (lower => more human)
double log_humanness_energy(const Chain* human_ref_seq, const char* seq, int n) {
    double energy = 0.0;
    for(int i=0; i<n; i++) {
        int idx = char_to_int(seq[i]);
        if(idx == -1) {
            fprintf(stderr, "Error: invalid character, %c, at position %d in sequence; penalizing as if it were 0.\n", seq[i], i);
            energy += ZERO_FREQ_PENALTY_LOG;
            continue;
        }
        double p = human_ref_seq->aas[i].elements[idx];
        if(p > EPSILON) energy -= log(p);
        else energy += ZERO_FREQ_PENALTY_LOG;
    }
    return(energy);
}

double linear_humanness_energy(const Chain* human_ref_seq, const char* seq, int n) {
    double energy = 0.0;
    for(int i=0; i<n; i++) {
        int idx = char_to_int(seq[i]);
        if(idx == -1) {
            fprintf(stderr, "Error: invalid character, %c, at position %d in sequence; penalizing as if it were 0.\n", seq[i], i);
            energy += ZERO_FREQ_PENALTY_LINEAR;
            continue;
        }
        double p = human_ref_seq->aas[i].elements[idx];
        energy -= p;
    }
    return(energy);
}

double property_distance_energy(const Chain* human_ref_seq, const char* seq, int n) {
    double t_energy = 0.0;
    for(int i=0; i<n; i++) {
        int idx = char_to_int(seq[i]);
        if(idx == -1) {
            fprintf(stderr, "Error: invalid character, %c, at position %d in sequence; penalizing as if it were 0.\n", seq[i], i);
            t_energy += ZERO_FREQ_PENALTY_PROPERTIES_DISTANCE;
            continue;
        }

        int* properties = PROPS_AA[idx];
        double sq_euclid_distance = 0.0;
        for(int j=0; j<N_PROPERTIES; j++) {
            double diff = (double)properties[j] - human_ref_seq->aas[i].properties[j];
            sq_euclid_distance += diff*diff;
        }
        t_energy += sq_euclid_distance;
    }
    return(t_energy);
}

double calculate_energy(const char* seq1, const char* seq2, int n) {
    return(hamming_distance(seq1, seq2, n));
}

/************* 
 * Next steps:
 * Adding more entropies
 * Adding multiple betas (infty for CDR's)
 * Maybe think of a way to ponder property_distance AND log prob?
 * 
 * 
 * 
 * 
 * 
 * 
 *************/

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