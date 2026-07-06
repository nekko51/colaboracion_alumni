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

/**
* Next steps:
* Adding more energies
* Adding multiple betas (infty for CDR's)
* Maybe think of a way to ponder property_distance AND log prob?
**/

/**
 * @brief Performs one Metropolis sweep over the sequence, attempting to mutate each position
 * 
 * @param murine_seq the murine sequence to be humanized
 * @param human_freqs the schrödinger human chain
 * @param beta 1/Temperature
 * @param n seq length
 * @param w_log weight for the log humanness energy
 * @param w_prop weight for the property distance energy
 */
void metropolis_sweep(char* murine_seq, const Chain* human_freqs, double beta, int n, double w_log, double w_prop) {
    for(int i=0; i<n; i++) {
        //A. Store the old State (the current amino acid at this position)

        //B. Propose a new State (a random, different amino acid)

        //C. Calculate delta_e for the Log Humanness Energy

        //D. Calculate delta_e for the Property Distance Energy

        //E. Combine and Weigh the Deltas to get the total_delta_e

        //F. The Metropolis Test: accept or reject the change, reverting back to original values if rejected

    }
}