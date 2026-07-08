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
        energy += 1 - p;
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

//calculates energies for a given sequence
Energy energy_calculation(const Chain* ref, const char* seq, int n) {
    Energy out;
    out.log_humanness = log_humanness_energy(ref, seq, n);
    out.property_distance = property_distance_energy(ref, seq, n);
    return(out);
}

void cleanup_metropolis(double* acceptance, Energy* energy_history, char** murine_history, int n_betas) {
    free(acceptance);
    free(energy_history);
    if (murine_history != NULL) {
        for (int i=0; i<n_betas+1; i++) {
            free(murine_history[i]);
        }
        free(murine_history);
    }
}

void cleanup_mega_metropolis(FILE* f, char** murine_seeds, int n_lines) {
    if (f != NULL) {
        fclose(f);
    }
    if (murine_seeds != NULL) {
        for (int i=0; i<n_lines; i++) {
            free(murine_seeds[i]);
        }
        free(murine_seeds);
    }
}

void print_metropolis_data_to_file(const char** seq_history, const Chain* reference, const Energy* energy_history, 
    const double* betas, int n_betas, const double* acceptance, int n_steps,
    double w_log, double w_prop, const char* filename) {
    //WIP
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
void metropolis_sweep(char* murine_seq, const Chain* human_ref_seq, double beta, double* acceptance, int n, double w_log, double w_prop) {
    int cnt = 0;
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
        if (total_delta_e < 0.0 || Random() < exp(-beta*total_delta_e)) {
            murine_seq[i] = new_AA;
            cnt++;
        }
    }
    *acceptance = (double)cnt/(double)n;
}

//returns 1 if there's an error, 0 if everything went smoothly -- char murine_seq[CHAINLEN+1];
int run_metropolis(char* murine_seq, const Chain* human_ref_seq, int n_steps, double* betas, int n_betas, char* filename) {
    //initial definitions & initializations
    double beta;
    double *acceptance = NULL;
    Energy* energy_history = NULL;
    char** murine_history = NULL;
    acceptance = malloc(n_steps*sizeof(double));
    if(acceptance == NULL) {fprintf(stderr, "Error: couldn't allocate memory for acceptance array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas); return(1);}

    energy_history = malloc( (n_betas+1)*sizeof(Energy) );
    if(energy_history == NULL) {fprintf(stderr, "Error: couldn't allocate memory for energy_history array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas);return(1);}

    murine_history = malloc( (n_betas+1)*sizeof(char*) );
    if(murine_history == NULL) {fprintf(stderr, "Error: couldn't allocate memory for murine_history *array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas);return(1);}
    for(int i=0; i<n_betas+1; i++) murine_history[i] = NULL;
    for(int i=0; i<n_betas+1; i++) {
        murine_history[i] = malloc( (CHAINLEN+1)*sizeof(char) );
        if(murine_history[i] == NULL) {fprintf(stderr, "Error: couldn't allocate memory for murine_history array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas);return(1);}
    }

    if(negative_chain(human_ref_seq) == 1) {fprintf(stderr, "Error: human_ref_Seq is invalid; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas);return(1);}

    //initial seed data
    strcpy(murine_history[0], murine_seq);
    energy_history[0] = energy_calculation(human_ref_seq, murine_history[0], CHAINLEN);

    //monte carlo sims
    for(int i=0; i<n_betas; i++) {
        beta = betas[i];
        for(int j=0; j<n_steps; j++) {
            metropolis_sweep(murine_seq, human_ref_seq, beta, &acceptance[j], CHAINLEN, WEIGHT_LOG, WEIGHT_PROP);
        }
        strcpy(murine_history[i+1], murine_seq);
        energy_history[i+1] = energy_calculation(human_ref_seq, murine_history[i+1], CHAINLEN);
    }

    for(int i=0; i<n_betas; i++) {
        print_metropolis_data_to_file(murine_history, human_ref_seq, energy_history, betas, n_betas, acceptance, n_steps, 
            WEIGHT_LOG, WEIGHT_PROP, filename);
    }

    cleanup_metropolis(acceptance, energy_history, murine_history, n_betas);
    return 0;
}

//runs metropolis using every single sequence in a given filename as seed n_metropolis times for each seed
void mega_metropolis(char* filename, int n_steps, double* betas, int n_betas, int n_metropolis) {
    /*initial declarations for cleanup*/
    FILE *f = NULL;
    char** murine_seeds = NULL;
    int n_lines = 0;
    char time_str[80];
    time_t now = time(NULL);

    /*open file and count number of lines*/
    f = fopen(filename, "r");
    if(f == NULL) {fprintf(stderr, "Error: couldn't open seed file %s; returning...\n", filename); cleanup_mega_metropolis(f, murine_seeds, n_lines); return;}
    char temp[CHAINLEN+1];
    while(1) {
        read_next_line(f, temp);
        if(temp[0] == '\0') break;
        n_lines++;
    }
    fclose(f);
    f = NULL;

    /*allocate necessary memory*/
    Chain human_ref_seq = file_megaAacids(SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES);
    murine_seeds = malloc(n_lines*sizeof(char*));
    if (murine_seeds == NULL) {fprintf(stderr, "Error: failed to allocate memory for *seed; returning...\n"); cleanup_mega_metropolis(f, murine_seeds, n_lines); return;}
    for(int i=0; i<n_lines; i++) murine_seeds[i] = NULL;
    for (int i=0; i<n_lines; i++) {
        murine_seeds[i] = malloc(CHAINLEN + 1);
        if (murine_seeds[i] == NULL) {fprintf(stderr, "Error: failed to allocate memory for seed %d; returning...\n", i); cleanup_mega_metropolis(f, murine_seeds, n_lines); return;}
    }

    /*read all seeds in file*/
    f = fopen(filename, "r");
    if(f == NULL) {fprintf(stderr, "Error: couldn't open seed file %s; returning...\n", filename); cleanup_mega_metropolis(f, murine_seeds, n_lines); return;}
    for(int i=0; i<n_lines; i++) {
        read_next_line(f, murine_seeds[i]);
    }
    fclose(f);
    f = NULL;
    
    /*create directories and run metropolis*/
    char metropolis_dir[MAX_STR_LEN];
    char batch_dir[MAX_STR_LEN];
    
    snprintf(metropolis_dir, sizeof(metropolis_dir), "%s%s", RESULTS, METROPOLIS);

    strftime(time_str, sizeof(time_str), "%Y-%m-%d_%H-%M-%S", localtime(&now));
    snprintf(batch_dir, sizeof(batch_dir), "%s%s_l%d_m%d_b%d", metropolis_dir, time_str, n_lines, n_metropolis, n_betas);
    MKDIR(batch_dir);

    char current_seed[CHAINLEN+1];
    for(int i=0; i<n_lines; i++) {
        printf("Running metropolis %d time(s) for seed %d/%d...\n", n_metropolis, i+1, n_lines);
        
        char seq_dir[MAX_STR_LEN];
        snprintf(seq_dir, sizeof(seq_dir), "%s/seq_%d", batch_dir, i+1);
        MKDIR(seq_dir);

        for(int j=0; j<n_metropolis; j++) {
            strcpy(current_seed, murine_seeds[i]);
            char filepath[MAX_STR_LEN];
            snprintf(filepath, sizeof(filepath), "%s/run_%d%s", seq_dir, j+1, TXT);
            run_metropolis(current_seed, &human_ref_seq, n_steps, betas, n_betas, filepath);
        }
    }

    /*
    hey, we're just two stray lines of completely functional code
    definitely code, see? we're started by a * so we're a pointer!
    */

    /*free memory*/
    cleanup_mega_metropolis(f, murine_seeds, n_lines);
}
