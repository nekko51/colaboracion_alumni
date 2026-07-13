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

//sq euclid distance between AA and schrödinger chain properties
double property_distance_from_ref(int aa_idx, const double* ref_props) {
    double sq_dist = 0.0;
    for(int j=0; j<N_PROPERTIES; j++) {
        double diff = (double)PROPS_AA[aa_idx][j] - ref_props[j];
        sq_dist += diff*diff;
    }
    return(sq_dist);
}

//sq euclid distance between seq and ref chain
double property_distance_energy(const Chain* human_ref_seq, const char* seq, int n) {
    double t_energy = 0.0;
    for(int i=0; i<n; i++) {
        int idx = char_to_int(seq[i]);
        if(idx == -1) {
            fprintf(stderr, "Error: invalid character, %c, at position %d in sequence; penalizing as if it were 0.\n", seq[i], i);
            t_energy += ZERO_FREQ_PENALTY_PROPERTIES_DISTANCE;
            continue;
        }
        t_energy += property_distance_from_ref(idx, human_ref_seq->aas[i].properties);
    }
    return(t_energy);
}

//sq euclid distance between aas in property space
double property_distance_between_aas(int aa_idx1, int aa_idx2) {
    double sq_dist = 0.0;
    for(int j=0; j<N_PROPERTIES; j++) {
        double diff = (double)PROPS_AA[aa_idx1][j] - (double)PROPS_AA[aa_idx2][j];
        sq_dist += diff*diff;
    }
    return(sq_dist);
}

double wanderer_penalty_energy(const char* seq, const char* original_seq, int n) {
    double t_energy = 0.0;
    for(int i=0; i<n; i++) {
        t_energy += property_distance_between_aas(char_to_int(seq[i]), char_to_int(original_seq[i]));
    }
    return(t_energy);
}

//calculates energies for a given sequence
Energy energy_calculation(const Chain* ref, const char* seq, const char* original_seq, int n) {
    Energy out;
    out.log_humanness = log_humanness_energy(ref, seq, n);
    out.property_distance = property_distance_energy(ref, seq, n);
    out.wanderer_penalty = wanderer_penalty_energy(seq, original_seq, n);
    return(out);
}

//calculates total energy given an energy and corresponding weights
double calculate_total_energy(const Energy* energies, double w_log, double w_prop, double w_penalty) {
    return(w_log*energies->log_humanness + w_prop*energies->property_distance + w_penalty*energies->wanderer_penalty);
}

    /**
    @brief generates betas matrix
    @param entropies double array of entropies for each chain position
    @param n_betas (rows) number of betas to be ran in metropolis
    @param n_entropies use CHAINLEN (columns) 
    @param scale_Factor \alpha in README.md
    @param epsilon use EPSILON (avoids division by 0)
    @param cooling_rate c_r in README.md -- SHOULD BE > 1
    @return betas matrix according to README.md -> Metropolis
    **/
void generate_betas(double** betas, int n_betas, double* entropies, int n_entropies, double scale_factor, double epsilon, double cooling_rate) {
    double k_i = 1.0;
    for(int i=0; i<n_betas; i++) {
        for(int j=0; j<n_entropies; j++) {
            betas[i][j] = k_i * scale_factor / (epsilon + entropies[j]);
        }
        k_i *= cooling_rate;
    }
}

double sigmoid(double x, double lambda) {
    return ( MAX_AMORTIG / (1.0+(MAX_AMORTIG - 1.0)*exp(-lambda * x)) );
}

/************************
* 
* Metropolis
*
************************/
void cleanup_metropolis(double* acceptance, Energy* energy_history, char** murine_history, int n_betas, double* per_beta_acceptance) {
    free(acceptance);
    free(per_beta_acceptance);
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
    double** betas, int n_betas, const double* acceptance, int n_sweeps,
    double w_log, double w_prop, double w_penalty, const char* filename) {
    FILE *f = fopen(filename, "w");
    if(f == NULL) {fprintf(stderr, "Error: couldn't open file %s; returning...\n", filename); return;}
    
    //separate file for betas matrix
    char beta_filename[MAX_STR_LEN];
    strncpy(beta_filename, filename, sizeof(beta_filename));
    beta_filename[sizeof(beta_filename)-1] = '\0';
    //find last dot position
    char *dot = strrchr(beta_filename, '.');
    if (dot != NULL) {
        *dot = '\0';
    }
    strncat(beta_filename, "_betas.txt", sizeof(beta_filename)-strlen(beta_filename)-1);

    FILE *beta_f = fopen(beta_filename, "w");
    if (beta_f != NULL) {
        fprintf(beta_f, "Beta matrix values (n_betas x CHAINLEN)\n\n");
        for (int i=0; i<n_betas; i++) {
            fprintf(beta_f, "Beta step %d:\n", i+1);
            for (int j=0; j<CHAINLEN; j++) {
                fprintf(beta_f, "%.*lf\t", DECIMAL_PRECISION, betas[i][j]);
            }
            fprintf(beta_f, "\n\n");
        }
        fclose(beta_f);
    }
    
    fprintf(f, "***************Metropolis data report version 1.0.2***************\n");

    fprintf(f, "\n\n");
    fprintf(f, "weight_log = %.*lf, \tweight_properties = %.*lf, \tweight_penalty = %.*lf\nn_betas = %d, \tn_sweeps = %d\n", 
        DECIMAL_PRECISION, w_log, DECIMAL_PRECISION, w_prop, DECIMAL_PRECISION, w_penalty, n_betas, n_sweeps);
    double mean, var;
    double med_acceptance[n_betas];
    for(int i=0; i<n_betas; i++) {
        med_acceptance[i] = acceptance[2*i];
    }
    med_var(med_acceptance, &mean, &var, n_betas);
    fprintf(f, "overall acceptance med: %.*lf, \toverall acceptance var: %.*lf\n", DECIMAL_PRECISION, mean, DECIMAL_PRECISION, var);//this is probably useless
    fprintf(f, "log humanness: %*.*lf, \tprop humanness: %*.*lf, \twanderer penalty: %*.*lf, \ttotal energy: %*.*lf\n", 
        ENERGY_PRECISION, DECIMAL_PRECISION, energy_history[0].log_humanness, ENERGY_PRECISION, DECIMAL_PRECISION, energy_history[0].property_distance,
        ENERGY_PRECISION, DECIMAL_PRECISION, energy_history[0].wanderer_penalty,
        ENERGY_PRECISION, DECIMAL_PRECISION, calculate_total_energy(&energy_history[0], w_log, w_prop, w_penalty));
    fprintf(f, "seed: %s", seq_history[0]);
    fprintf(f, "\n\n");
    // fprtinf(f, "Schrödinger reference chain: ")
    for(int i=0; i<n_betas; i++) {
        double min_beta, max_beta, avg_beta = 0.0;
        minmax(betas[i], &max_beta, &min_beta, CHAINLEN);
        for(int j=0; j<CHAINLEN; j++) {
            avg_beta += betas[i][j];
        }
        avg_beta /= CHAINLEN;
        fprintf(f, "***** beta (min/avg/max) = %*.*lf / %*.*lf / %*.*lf (%d/%d) \t-\t acceptance med: %.*lf%%, \tacceptance var: %.*lf*****\n",
            BETA_PRECISION, DECIMAL_PRECISION, min_beta,
            BETA_PRECISION, DECIMAL_PRECISION, avg_beta,
            BETA_PRECISION, DECIMAL_PRECISION, max_beta,
            i+1, n_betas, ACCEPTANCE_PRECISION, 100.0*acceptance[2*i], ACCEPTANCE_PRECISION, acceptance[2*i+1]);
        fprintf(f, "delta_e: %*.*lf, \tlog humanness: %*.*lf, \tprop humanness: %*.*lf, \twanderer penalty: %*.*lf, \ttotal energy: %*.*lf\n", 
            ENERGY_PRECISION, DECIMAL_PRECISION, calculate_total_energy(&energy_history[i+1], w_log, w_prop, w_penalty) - calculate_total_energy(&energy_history[i], w_log, w_prop, w_penalty), 
            ENERGY_PRECISION, DECIMAL_PRECISION, energy_history[i+1].log_humanness, ENERGY_PRECISION, DECIMAL_PRECISION, energy_history[i+1].property_distance,
            ENERGY_PRECISION, DECIMAL_PRECISION, energy_history[i+1].wanderer_penalty,
            ENERGY_PRECISION, DECIMAL_PRECISION, calculate_total_energy(&energy_history[i+1], w_log, w_prop, w_penalty));
        fprintf(f, "resulting chain: %s\n", seq_history[i+1]);
    }

    fprintf(f, "\n\n***************Metropolis data report end***************");
    fclose(f);
}


void metropolis_sweep(char* murine_seq, const char* original_murine_seq, const Chain* human_ref_seq, double* local_beta, 
    double* acceptance, int chainlen, double w_log, double w_prop, double w_penalty) {
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
    int cnt = 0;
    for(int i=0; i<chainlen; i++) {
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
        Energy delta_e;
        double old_p = human_ref_seq->aas[i].elements[old_AA_idx];
        double old_log_contrib;
        if(old_p > EPSILON) old_log_contrib = -log(old_p);
        else old_log_contrib = ZERO_FREQ_PENALTY_LOG;

        double new_p = human_ref_seq->aas[i].elements[new_AA_idx];
        double new_log_contrib;
        if(new_p > EPSILON) new_log_contrib = -log(new_p);
        else new_log_contrib = ZERO_FREQ_PENALTY_LOG;

        delta_e.log_humanness = new_log_contrib - old_log_contrib;

        //D. Calculate delta_e for the Property Distance Energy
        double old_prop_dist = property_distance_from_ref(old_AA_idx, human_ref_seq->aas[i].properties);
        double new_prop_dist = property_distance_from_ref(new_AA_idx, human_ref_seq->aas[i].properties);
        delta_e.property_distance = new_prop_dist - old_prop_dist;

        //E. Calculate delta_e for penalty term
        int original_AA_idx = char_to_int(original_murine_seq[i]);
        /*3 possible cases:
        * if both are different from original, penalties cancel each other
        * if new isn't original and old one is, a penalty is added
        * if old isn't original and new one is, penalty is reverted
        */
        double old_penalty_dist = property_distance_between_aas(old_AA_idx, original_AA_idx);
        if(old_AA_idx != original_AA_idx) {
            old_penalty_dist += AA_MUTATION_PENALTY;
        }
        double new_penalty_dist = property_distance_between_aas(new_AA_idx, original_AA_idx);
        if(new_AA_idx != original_AA_idx) {
            new_penalty_dist += AA_MUTATION_PENALTY;
        }

        delta_e.wanderer_penalty = new_penalty_dist - old_penalty_dist;
        
        //F. Combine and Weigh the Deltas to get the total_delta_e
        double total_delta_e = calculate_total_energy(&delta_e, w_log, w_prop, w_penalty);

        //G. The Metropolis Test: accept or reject the change, reverting back to original values if rejected
        if (total_delta_e < 0.0 || Random() < exp(-local_beta[i]*total_delta_e)) {
            murine_seq[i] = new_AA;
            cnt++;
        }
    }
    *acceptance = (double)cnt/(double)chainlen;
}

//returns 1 if there's an error, 0 if everything went smoothly -- runs n_sweeps metropolis sweeps
int run_metropolis(char* murine_seq, const Chain* human_ref_seq, int n_sweeps, double** betas, int n_betas, char* filename) {
    //initial definitions & initializations
    double local_beta[CHAINLEN];
    double *acceptance = NULL;
    double *per_beta_acceptance = NULL;
    char original_murine_seq[CHAINLEN+1];
    strcpy(original_murine_seq, murine_seq);
    Energy* energy_history = NULL;
    char** murine_history = NULL;
    acceptance = malloc(n_sweeps*sizeof(double));
    if(acceptance == NULL) {fprintf(stderr, "Error: couldn't allocate memory for per-sweep acceptance array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas, per_beta_acceptance); return(1);}
    per_beta_acceptance = malloc(2*n_betas*sizeof(double));
    if(per_beta_acceptance == NULL) {fprintf(stderr, "Error: couldn't allocate memory for per-beta acceptance array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas, per_beta_acceptance); return(1);}

    energy_history = malloc( (n_betas+1)*sizeof(Energy) );
    if(energy_history == NULL) {fprintf(stderr, "Error: couldn't allocate memory for energy_history array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas, per_beta_acceptance);return(1);}

    murine_history = malloc( (n_betas+1)*sizeof(char*) );
    if(murine_history == NULL) {fprintf(stderr, "Error: couldn't allocate memory for murine_history *array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas, per_beta_acceptance);return(1);}
    for(int i=0; i<n_betas+1; i++) murine_history[i] = NULL;
    for(int i=0; i<n_betas+1; i++) {
        murine_history[i] = malloc( (CHAINLEN+1)*sizeof(char) );
        if(murine_history[i] == NULL) {fprintf(stderr, "Error: couldn't allocate memory for murine_history array; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas, per_beta_acceptance);return(1);}
    }

    if(negative_chain(human_ref_seq) == 1) {fprintf(stderr, "Error: human_ref_Seq is invalid; returning 1 in metropolis..."); cleanup_metropolis(acceptance, energy_history, murine_history, n_betas, per_beta_acceptance);return(1);}

    //initial seed data
    strcpy(murine_history[0], murine_seq);
    energy_history[0] = energy_calculation(human_ref_seq, murine_history[0], original_murine_seq, CHAINLEN);

    //monte carlo sims
    for(int i=0; i<n_betas; i++) {
        for(int j=0; j<CHAINLEN; j++) local_beta[j] = betas[i][j];

        double progress = (double)(i+1) / (double)n_betas;
        double current_target = STARTING_TARGET_ACCEPTANCE*(1.0-progress);

        int local_success[CHAINLEN];

        double last_acceptance = current_target;//we begin by assuming we're in our target acceptance

        for(int j=0; j<n_sweeps; j++) {
            double deviation = last_acceptance - current_target;
            double factor = sigmoid(deviation, LAMBDA);

            double amortig_beta[CHAINLEN];
            for(int k=0; k<CHAINLEN; k++) amortig_beta[k]=factor*local_beta[k];

            metropolis_sweep(murine_seq, original_murine_seq, human_ref_seq, amortig_beta, &acceptance[j], CHAINLEN, WEIGHT_LOG, WEIGHT_PROP, WEIGHT_PENALTY);
            last_acceptance = acceptance[j];
        }
        strcpy(murine_history[i+1], murine_seq);
        energy_history[i+1] = energy_calculation(human_ref_seq, murine_history[i+1], original_murine_seq, CHAINLEN);
        med_var(acceptance, &per_beta_acceptance[2*i], &per_beta_acceptance[2*i+1], n_sweeps);
    }

    print_metropolis_data_to_file((const char**)murine_history, human_ref_seq, energy_history, betas, n_betas, 
        per_beta_acceptance, n_sweeps, WEIGHT_LOG, WEIGHT_PROP, WEIGHT_PENALTY, filename);
    
    cleanup_metropolis(acceptance, energy_history, murine_history, n_betas, per_beta_acceptance);
    return 0;
}

//runs metropolis using every single sequence in a given filename as seed n_metropolis times for each seed
void mega_metropolis(char* murine_seeds_filename, char* human_filename, int n_human_lines, int n_sweeps, double** betas, int n_betas, int n_metropolis) {
    /*initial declarations for cleanup*/
    FILE *f = NULL;
    char** murine_seeds = NULL;
    int n_lines = 0;
    char time_str[80];
    time_t now = time(NULL);

    /*open file and count number of lines*/
    f = fopen(murine_seeds_filename, "r");
    if(f == NULL) {fprintf(stderr, "Error: couldn't open seed file %s; returning...\n", murine_seeds_filename); cleanup_mega_metropolis(f, murine_seeds, n_lines); return;}
    char temp[CHAINLEN+1];
    while(1) {
        read_next_line(f, temp);
        if(temp[0] == '\0') break;
        n_lines++;
    }
    fclose(f);
    f = NULL;

    /*allocate necessary memory*/
    murine_seeds = malloc(n_lines*sizeof(char*));
    if (murine_seeds == NULL) {fprintf(stderr, "Error: failed to allocate memory for *seed; returning...\n"); cleanup_mega_metropolis(f, murine_seeds, n_lines); return;}
    for(int i=0; i<n_lines; i++) murine_seeds[i] = NULL;
    for (int i=0; i<n_lines; i++) {
        murine_seeds[i] = malloc(CHAINLEN + 1);
        if (murine_seeds[i] == NULL) {fprintf(stderr, "Error: failed to allocate memory for seed %d; returning...\n", i); cleanup_mega_metropolis(f, murine_seeds, n_lines); return;}
    }

    /*read all seeds in file*/
    f = fopen(murine_seeds_filename, "r");
    if(f == NULL) {fprintf(stderr, "Error: couldn't open seed file %s; returning...\n", murine_seeds_filename); cleanup_mega_metropolis(f, murine_seeds, n_lines); return;}
    for(int i=0; i<n_lines; i++) {
        read_next_line(f, murine_seeds[i]);
    }
    fclose(f);
    f = NULL;
    
    /*create directories and run metropolis*/
    char metropolis_dir[MAX_STR_LEN];
    char batch_dir[MAX_STR_LEN];
    
    snprintf(metropolis_dir, sizeof(metropolis_dir), "%s%s", RESULTS, METROPOLIS);

    strftime(time_str, sizeof(time_str), "%Y-%m-%d_%H_%M_%S", localtime(&now));
    snprintf(batch_dir, sizeof(batch_dir), "%s%s_l%d_m%d_b%d", metropolis_dir, time_str, n_lines, n_metropolis, n_betas);
    MKDIR(batch_dir);

    Chain human_ref_seq = file_megaAacids(human_filename, n_human_lines);

    int workers = (int) ((MAX_WORKER_PERCENTAGE/100.0) * omp_get_num_procs());
    if(workers < 1) workers = 1;
    int global_error = 0;
    #pragma omp parallel for num_threads(workers) schedule(dynamic)
    for(int i=0; i<n_lines; i++) {
        //openmp doesn't allow breaks, and the whole for has to be ran; if there's an error, threads just skip everything and end ASAP
        if(global_error) continue;

        printf("Running metropolis %d time(s) for seed %d/%d...\n", n_metropolis, i+1, n_lines);
        char current_seed[CHAINLEN+1];
        char seq_dir[MAX_STR_LEN];
        char filepath[MAX_STR_LEN];
        snprintf(seq_dir, sizeof(seq_dir), "%s/seq_%d", batch_dir, i+1);
        MKDIR(seq_dir);

        int okay = 0;
        int j;
        for(j=0; j<n_metropolis && okay == 0; j++) {
            strcpy(current_seed, murine_seeds[i]);
            snprintf(filepath, sizeof(filepath), "%s/run_%d%s", seq_dir, j+1, TXT);

            okay = run_metropolis(current_seed, &human_ref_seq, n_sweeps, betas, n_betas, filepath);
        }
        if(okay != 0) {
            fprintf(stderr, "Error: run_metropolis failed in interation %d; returning...\n", j);
            global_error = 1;
        }
    }

    //safe error handling, oustide parallel for loop
    if(global_error) {
        cleanup_mega_metropolis(f, murine_seeds, n_lines); 
        return;
    }

    /*
    hey, we're just two stray lines of completely functional code
    definitely code, see? we're started by a * so we're a pointer!
    */

    /*free memory*/
    cleanup_mega_metropolis(f, murine_seeds, n_lines);
}