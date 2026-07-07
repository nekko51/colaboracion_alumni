#include "head.h"

/************************
* 
* Energies
*
************************/
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

//UNFINISHED (WIP)!
//n = 0 for chain of zeroes, n = 1 for schrödinger mouse chain, any other n will result in error, -1.0 chain
Chain generate_murine_seed(seed n) {
    Chain out;
    switch(n) {
        case zeroes:
            memset(&out, 0, sizeof(Chain));
            break;
        case megachain:
            out = file_megaAacids(SEQS FILE_L_MOUSE TXT, L_MOUSE_N_LINES);
            break;
        case error:
        default:
            for(int i=0; i<CHAINLEN; i++) {
                for(int j=0; j<N_AACIDS; j++) out.aas[i].elements[j] = -1.;
                for(int j=0; j<N_PROPERTIES; j++) out.aas[i].properties[j] = -1.;
            }
            break;
    }
    return(out);
}

/**
 * @brief Performs one Metropolis sweep over the sequence, attempting to mutate each position
 * 
 * @param murine_seq the murine sequence to be humanized
 * @param human_freqs the schrödinger human chain
 * @param beta 1/Temperature
 * @param n seq length
 * @param w_log weight for the log humanness energy
 * @param w_prop weight for the property distance energy
 * 
 * Next steps:
 * Adding more energies/changing weight system
 * Adding multiple betas (infty for CDR's)
 */
void metropolis_sweep(char* murine_seq, const Chain* human_ref_seq, double beta, int n, double w_log, double w_prop) {
    for(int i=0; i<n; i++) {
        //A. Store the old State (the current amino acid at this position)
        char old_AA = murine_seq[i];
        int old_AA_idx = char_to_int(old_AA);
        if (old_AA_idx == -1) {
            fprintf(stderr, "Error: AA sampled isn't a valid AA, '%c' (ASCII %d); continuing loop\n", old_AA, (int) old_AA);
            continue;
        }

        //B. Propose a new State (a random, different amino acid)
        int new_AA_idx;
        do {
            new_AA_idx = (int) (N_AACIDS*Random());
        } while(new_AA_idx == old_AA_idx);
        char new_AA = AMINOACIDS[new_AA_idx];

        //C. Calculate delta_e for the Log Humanness Energy
        double old_p = human_ref_seq->aas[i].elements[old_AA_idx];
        double old_log_contrib;
        if(old_p > EPSILON) old_log_contrib = -log(old_p);
        else old_log_contrib = ZERO_FREQ_PENALTY_LOG;

        double new_p = human_ref_seq->aas[i].elements[new_AA_idx];
        double new_log_contrib;
        if(new_p > EPSILON) new_log_contrib = -log(new_p);
        else new_log_contrib = ZERO_FREQ_PENALTY_LOG;

        double delta_log = new_log_contrib - old_log_contrib;

        //D. Calculate delta_e for the Property Distance Energy
        double old_prop_dist = 0.0, new_prop_dist = 0.0;
        for(int j=0; j<N_PROPERTIES; j++) {
            double old_diff = (double)PROPS_AA[old_AA_idx][j] - human_ref_seq->aas[i].properties[j];
            double new_diff = (double)PROPS_AA[new_AA_idx][j] - human_ref_seq->aas[i].properties[j];
            old_prop_dist += old_diff*old_diff;
            new_prop_dist += new_diff*new_diff;
        }

        double delta_prop = new_prop_dist - old_prop_dist;

        //E. Combine and Weigh the Deltas to get the total_delta_e
        double total_delta_e = w_log*delta_log + w_prop*delta_prop;

        //F. The Metropolis Test: accept or reject the change, reverting back to original values if rejected
        if (total_delta_e < 0.0 || Random() < exp(-beta*total_delta_e)) murine_seq[i] = new_AA;
    }
}