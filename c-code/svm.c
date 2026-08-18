#include "head.h"



// @return exp( - SVD_DISTANCE_PROP_FACTOR * chain_sq_distance(a, b))
double gaussian_radial_basis_function(Chain a, Chain b) {
    return exp( - SVM_KERNEL_PROP_FACTOR * chain_sq_distance(a, b) );
}


double sigmoid_function(Chain a, Chain b) {
    return tanh( SVM_KERNEL_PROP_FACTOR * (chain_dot_product(a, b) + SVM_KERNEL_LAG));
}

double polynomial_inhom_function(Chain a, Chain b) {
    double lin_poly = chain_dot_product(a, b) + SVM_KERNEL_LAG;
    double ret = 1.;
    for (int i = 0; i < SVM_KERNEL_POW; i++) {
        ret *= lin_poly;
    }
    return ret;
}

/**
 * 's' for sigmoid (tanh)
 * 'p' for polynomial
 * any other one will be for gaussian (RBF)
 */
double kernel(Chain a, Chain b, char kernel_char) {
    switch (kernel_char) {
    case 's':   return sigmoid_function(a, b);
    case 'p':   return polynomial_inhom_function(a, b);
    default:    return gaussian_radial_basis_function(a, b);
    };
}

// // @brief [a_i], [b_j], takes the index-th element in the concatenation of *c1 and *c2
// // @returns -1 if it was taken from the 1st list, 1 if from the 2nd
// int get_chain_from_yuxtaposed_arrays(Chain *c1, int n_c1, Chain *c2, int n_c2, int index, Chain* out) {
//     if (index >= n_c1) {
//         *out = c2[index - n_c1];
//         return 1; // return integer if taken from the 2nd list
//     }
//     *out = c1[index];
//     return -1;
// }


/**
 * @brief computes (signed) gram matrix for a series of learn chains
 * @returns a matrix of size n_chains * n_chains
 * 
 * @param kernel_char choose the kernel for the outut matrix
 * 
 * - 's' for sigmoid (tanh): tanh( k * (a·b + r) )
 * 
 * - 'p' for polynomial: ( a·b + r )^d
 * 
 * - any other (standard to g) for gaussian RBF: exp( - k * ||a - b|| )
 */
void svm_gram_matrix(Chain *chains, int n_chains, double **out, char kernel_char) {
    printf("starting gram matrix computation...\n");
    for (int i = 0; i < n_chains; i++) {
        for (int j = 0; j < n_chains; j++) {
            out[i][j] = kernel(chains[i], chains[j], kernel_char);

        }
        printf("progress: %5.2f%%\r", (double)(i+1)/(double)n_chains * 100.);
        fflush(stdout);
    }
    printf("finished gram matrix!\n");
}

/**
 * decision function
 * f(x) = b + \sum_i \lambda_i s_i kern(x_i, x)
 * 
 * bias
 * b = s_p - \sum_i \lambda_i s_i kern(x_i, x_p) 
 * where \lambda_p > 0
 * 
 */

double compute_bias(double *lambda, int *signs, Chain *chains, int n_chains, char kernel_char) {
    int idx = -1;
    for (int i = 0; i < n_chains; i++) {
        if (lambda[i] > EPSILON) {
            idx = i;
            break;
        }
    }
    if (idx == -1) {
        printf("error trying to find support vector for bias: no support vectors were found");
        return 0.;
    }

    double bias = signs[idx];
    for (int i = 0; i < n_chains; i++) {
        bias -= lambda[i] * signs[i] * kernel(chains[i], chains[idx], kernel_char);
    }
    return bias;
}

double decision_function(Chain x_input, double bias, double *lambda, int *signs, Chain *chains, int n_chains, char kernel_char) {
    double out = bias;
    for (int i = 0; i < n_chains; i++) {
        out += lambda[i] * signs[i] * kernel(x_input, chains[i], kernel_char);
    }
    return out;
}


/*********************************************************************************************************************************/
/*******************************************************simulated annealing*******************************************************/
/*********************************************************************************************************************************/



SvmMove svm_propose_move(double *v, int v_dim, double eps, int *signs, double C_limit) {
    SvmMove move;
    int intentos = 0;
    double min_delta, max_delta;

    do {
        move.i = (int)(Random() * v_dim);
        if (move.i >= v_dim) move.i = v_dim - 1;
        do {
            move.j = (int)(Random() * v_dim);
            if (move.j >= v_dim) move.j = v_dim - 1;
        } while (move.i == move.j);

        double min_di = -v[move.i];
        double max_di = C_limit - v[move.i];
        double min_dj, max_dj;

        if (signs[move.i] == signs[move.j]) {
            min_dj = v[move.j] - C_limit;
            max_dj = v[move.j];
        } else {
            min_dj = -v[move.j];
            max_dj = C_limit - v[move.j];
        }

        min_delta = min_di > min_dj ? min_di : min_dj;
        max_delta = max_di < max_dj ? max_di : max_dj;
        
        double eps_min = -eps > min_delta ? -eps : min_delta;
        double eps_max = eps < max_delta ? eps : max_delta;

        if (eps_max > eps_min) {
            move.delta = eps_min + Random() * (eps_max - eps_min);
            return move;
        }
        
        intentos++;
        if (intentos > 1000) {
            move.delta = 0.0;
            return move;
        }
    } while (1); 
}


/**
 * maximize 
 * \sum_i \lambda_i - \frac{1}{2} \sum_i \sum_j \lambda_i \lambda_j s_i s_j * kern(x_i, x_j)
 * 
 * with
 * \lambda_i \ge 0
 * \sum_i \lambda_i s_i = 0
 */

double svm_energy(double **gram_matrix, int *signs, double *lambda, int dim) {
    double sum = 0.;
    for (int i = 0; i < dim; i++) {
        sum -= lambda[i];
        for (int j = 0; j < dim; j++) {
            sum += gram_matrix[i][j] * lambda[i] * lambda[j] * signs[i] * signs[j] / 2;
        }
    }
    return sum;
}

void svm_change(double *v, double *v_, int v_dim, double eps, int *signs, double C_limit) {
    for (int k = 0; k < v_dim; k++) v_[k] = v[k];

    int i, j;
    double min_delta, max_delta, delta;

    do {
        i = (int)(Random() * v_dim);
        if (i >= v_dim) i = v_dim - 1;
        do {
            j = (int)(Random() * v_dim);
            if (j >= v_dim) j = v_dim - 1;
        } while (i == j);

        double min_di = -v[i];
        double max_di = C_limit - v[i];
        double min_dj, max_dj;

        if (signs[i] == signs[j]) {
            min_dj = v[j] - C_limit;
            max_dj = v[j];
        } else {
            min_dj = -v[j];
            max_dj = C_limit - v[j];
        }

        min_delta = min_di > min_dj ? min_di : min_dj;
        max_delta = max_di < max_dj ? max_di : max_dj;
        
        double eps_min = -eps > min_delta ? -eps : min_delta;
        double eps_max = eps < max_delta ? eps : max_delta;

        if (eps_max > eps_min) {
            delta = eps_min + Random() * (eps_max - eps_min);
            v_[i] = v[i] + delta;
            v_[j] = v[j] - delta * signs[i] * signs[j];
            break;
        }
    } while (1); 
}

int svm_metropolis(double e_old, double e_new, double beta) {
    double ce = exp(-beta*(e_new - e_old));
    if (ce > Random()) {
        return 1;
    }
    return 0;
}

double svm_initial_beta(double *initial_lambdas, double eps, int iterations, double **gram_matrix, int *signs, int dim, double C_limit) {
    printf("starting initial beta computation...\n");
    
    double *F = calloc(dim, sizeof(double));
    for (int k = 0; k < dim; k++) {
        for (int m = 0; m < dim; m++) {
            F[k] += initial_lambdas[m] * signs[m] * gram_matrix[k][m];
        }
    }

    double sum = 0.0;
    
    for (int it = 0; it < iterations; it++) {
        SvmMove move = svm_propose_move(initial_lambdas, dim, eps, signs, C_limit);
        
        if (move.delta != 0.0) {
            double d = move.delta;
            int i = move.i;
            int j = move.j;
            
            double delta_E = -d * (1.0 - signs[i] * signs[j]) 
                             + d * signs[i] * (F[i] - F[j]) 
                             + 0.5 * d * d * (gram_matrix[i][i] + gram_matrix[j][j] - 2.0 * gram_matrix[i][j]);
            
            sum += fabs(delta_E);
        }
        
        if (it % 100 == 0 || it == iterations - 1) {
            printf("progress:\t%5.2f%%\r", (double)(it+1)/(double)iterations * 100.);
            fflush(stdout);
        }
    }
    
    free(F);
    
    double ret = 1.0;
    if (sum > EPSILON) {
        ret = (double)iterations / sum;
    }
    
    printf("initial beta is b = %.3g\n\n", ret);
    return ret;
}

void svm_initialize_lambdas(double *lambdas, int *signs, int dim, double C_limit) {
    double sum_pos = 0.0;
    double sum_neg = 0.0;

    for (int i = 0; i < dim; i++) {
        lambdas[i] = Random() * (C_limit * 0.1);
        if (signs[i] == 1) sum_pos += lambdas[i];
        else               sum_neg += lambdas[i];
    }

    if (sum_pos > sum_neg && sum_pos > 0.0) {
        for (int i = 0; i < dim; i++) {
            if (signs[i] == 1) lambdas[i] *= (sum_neg / sum_pos);
        }
    } else if (sum_neg > 0.0) {
        for (int i = 0; i < dim; i++) {
            if (signs[i] == -1) lambdas[i] *= (sum_pos / sum_neg);
        }
    }
    printf("lambdas initialized!\n\n");
}

/**
 * @brief computes simulated annealing for lagrange multipliers to optimize dual larangian problem for SVM
 * @param gram_matrix kernel matrix for each data point
 * @param signs tags of the data (-1 for murine, +1 for human)
 * @param dimension number of lambdas = number of data points
 * @param epsilon maximum distance of random leaps
 * @param betas beta for each iteration of betas
 * @param itertions_per_beta number of metropolis steps to the vector for each beta in list
 * @param print_to_file 0 for no printing; 1 for printing log of energies, lambda vectors and acceptance
 * @return parameter ``lambdas`` is used as first lambda vector used and output lambda after annealing
 * more to do: acceptance, outputting energy, printing to file
 */
void svm_annealing(double **gram_matrix, int *signs, int dimension, double epsilon, double *lambdas, double *betas, int n_betas, int iterations_per_beta, int print_to_file, double C_limit, int run_id) {
    printf("starting simulated annealing...\n");
    
    double energy_old = svm_energy(gram_matrix, signs, lambdas, dimension);
    int acceptance_count;
    
    double *F = calloc(dimension, sizeof(double));
    for (int k = 0; k < dimension; k++) {
        for (int m = 0; m < dimension; m++) {
            F[k] += lambdas[m] * signs[m] * gram_matrix[k][m];
        }
    }

    FILE *sum_file = NULL;
    FILE *log_file = NULL;
    char time_str[80];
    time_t now = time(NULL);
    strftime(time_str, sizeof(time_str), "%Y-%m-%d_%H:%M:%S", localtime(&now));
    char run_dir[MAX_STR_LEN];
    sprintf(run_dir, "results/svm/annealing/batches/betas%05d_iterations%05d_c-limit%g/%s-%03d", n_betas, iterations_per_beta, C_limit, time_str, run_id);
    
    if (print_to_file) {
    
        mkdir_p(run_dir);

        char filename[3*MAX_STR_LEN];
        char summry_file_name[3*MAX_STR_LEN];
        char log_file_name[3*MAX_STR_LEN];
        sprintf(filename, "%s/test-results.txt", run_dir);
        sprintf(summry_file_name, "%s/summary.txt", run_dir);
        sprintf(log_file_name, "%s/full-log.txt", run_dir);

        sum_file = get_file(summry_file_name, "w");
        log_file = get_file(log_file_name, "w");

        fprintf(sum_file, "B-idx\tBeta\tEnergy\tAcceptance\n");
        fprintf(log_file, "Line\tB-idx\tBeta\tIt-idx\tEnergy\tAcceptance~\n");
    }

    printf("progress:\t%5.2f%%\r", 0.);
    fflush(stdout);



    for (int beta_idx = 0; beta_idx < n_betas; beta_idx++) {
        acceptance_count = 0;

        for (int iteration_idx = 0; iteration_idx < iterations_per_beta; iteration_idx++) {
            
            SvmMove move = svm_propose_move(lambdas, dimension, epsilon, signs, C_limit);

            if (move.delta != 0.0) {
                double d = move.delta;
                int i = move.i;
                int j = move.j;
                
                double delta_E = -d * (1.0 - signs[i] * signs[j]) 
                                 + d * signs[i] * (F[i] - F[j]) 
                                 + 0.5 * d * d * (gram_matrix[i][i] + gram_matrix[j][j] - 2.0 * gram_matrix[i][j]);

                if (delta_E < 0.0 || exp(-betas[beta_idx] * delta_E) > Random()) {
                    
                    lambdas[i] += d;
                    lambdas[j] -= d * signs[i] * signs[j];
                    
                    if (lambdas[i] < 0.0) lambdas[i] = 0.0;
                    if (lambdas[i] > C_limit) lambdas[i] = C_limit;
                    if (lambdas[j] < 0.0) lambdas[j] = 0.0;
                    if (lambdas[j] > C_limit) lambdas[j] = C_limit;

                    energy_old += delta_E;
                    acceptance_count++;

                    for (int k = 0; k < dimension; k++) {
                        F[k] += d * signs[i] * (gram_matrix[k][i] - gram_matrix[k][j]);
                    }
                }
            }

            if ((iteration_idx + 1) % 100 == 0 || iteration_idx == iterations_per_beta - 1) {
                printf("progress:\t%5.2f%%\r", (double)(beta_idx * iterations_per_beta + iteration_idx + 1) / (double)(n_betas * iterations_per_beta) * 100.);
                fflush(stdout);
            }

            if (print_to_file) {
                fprintf(log_file, "%d\t%d\t%lf\t%d\t%lf\t%g\n", 
                        iteration_idx + beta_idx * iterations_per_beta, 
                        beta_idx, betas[beta_idx], iteration_idx, 
                        energy_old, (double)acceptance_count / (iteration_idx + 1));
            }
        }

        energy_old = svm_energy(gram_matrix, signs, lambdas, dimension);
        for (int k = 0; k < dimension; k++) {
            F[k] = 0.0;
            for (int m = 0; m < dimension; m++) {
                F[k] += lambdas[m] * signs[m] * gram_matrix[k][m];
            }
        }

        if (print_to_file) {
            fprintf(sum_file, "%d\t%lf\t%lf\t%lf\n", 
                    beta_idx + 1, betas[beta_idx], energy_old, 
                    (double)acceptance_count / (double)iterations_per_beta);
        }

        if ((double)acceptance_count / (double)iterations_per_beta > 0.5) {
            epsilon *= 1.05; 
        } else {
            epsilon *= 0.95;
        }
    }

    if (print_to_file) {
        char res_file_name[3*MAX_STR_LEN];
        sprintf(res_file_name, "%s/result-lambdas.txt", run_dir);
        FILE *res_file = get_file(res_file_name, "w");

        fprintf(res_file, "%g\t", energy_old);
        for(int i = 0; i < dimension; i++) {
            fprintf(res_file, "%g\t", lambdas[i]);
        }

        fclose(res_file);
        fclose(sum_file);
        fclose(log_file);
    }

    free(F);
    printf("finished simulated annealing!\n\n");
}

/**
 * @brief Reads an aggregated CSV file, finds the row with the minimum energy,
 *        and loads the corresponding lambda values.
 *
 * @param filename The path to the aggregated CSV file.
 * @param best_lambdas A pointer to a double array (size n_data_points) to store the resulting lambdas.
 * @param n_data_points The number of data points (and lambdas) to read.
 * @return 0 on success, -1 on file error, -2 if no data found.
 */
void get_best_lambdas_from_csv(double* best_lambdas, int n_data_points) {
    const char* filename = "results/svm/annealing/aggregated_svm_lambdas.csv";
    FILE *f = get_file((char*)filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Error: Could not open aggregated lambda file %s\n", filename);
        return;
    }

    char line[MAX_STR_LEN * 100]; // Buffer to hold a full line from the CSV
    char best_line[MAX_STR_LEN * 100] = "";
    double min_energy = 1e30; // Initialize with a very large number

    // Skip header line
    if (fgets(line, sizeof(line), f) == NULL) {
        fprintf(stderr, "Error: CSV file %s is empty or unreadable.\n", filename);
        fclose(f);
        return;
    }

    // Find the line with the minimum energy
    while (fgets(line, sizeof(line), f) != NULL) {
        double current_energy;
        // sscanf to extract the third value (energy)
        if (sscanf(line, "%*[^,],%*[^,],%lf", &current_energy) == 1) {
            if (current_energy < min_energy) {
                min_energy = current_energy;
                strcpy(best_line, line);
            }
        }
    }
    fclose(f);

    if (strlen(best_line) == 0) {
        fprintf(stderr, "No data found in %s\n", filename);
        return;
    }

    printf("Found best run with energy: %g\n", min_energy);

    // Parse the best line to extract lambdas
    char *token = strtok(best_line, ",");
    // Skip n_betas, n_iterations, and energy
    for (int i = 0; i < 3; ++i) {
        if (token == NULL) {
            fprintf(stderr, "Error: Malformed CSV line.\n");
            return;
        }
        token = strtok(NULL, ",");
    }

    // Read the lambda values
    for (int i = 0; i < n_data_points; ++i) {
        if (token != NULL) {
            best_lambdas[i] = atof(token);
            token = strtok(NULL, ",");
        } else {
            fprintf(stderr, "Error: Not enough lambda values in CSV line for n_data_points=%d.\n", n_data_points);
            break;
        }
    }
    printf("Successfully loaded best lambda vector.\n");
}

/********************************************************************************************************************************/
/***************************************************end of simulated annealing***************************************************/
/********************************************************************************************************************************/

Int2 test_results_to_file(int learn_dims, double *lambdas, int *tags, Chain *learn_chs, Chain *test_chs, int n_test_chs, char kernel_char, int correct_sign) {
    printf("starting tests...\n");
    double bias = compute_bias(lambdas, tags, learn_chs, learn_dims, kernel_char);

    char time_str[80];
    time_t now = time(NULL);
    strftime(time_str, sizeof(time_str), "%Y-%m-%d_%H:%M:%S", localtime(&now));
    char run_dir[MAX_STR_LEN];
    sprintf(run_dir, "results/svm/decision-function-tests/kernel-%c/%s", kernel_char, time_str);
    mkdir_p(run_dir);
    printf("directory '%s' created!\n", time_str);

    char filename[2*MAX_STR_LEN];
    sprintf(filename, "%s/test-results.txt", run_dir);

    FILE *f = get_file(filename, "w");
    double decision_f;
    int correct = 0;

    for (int i = 0; i < n_test_chs; i++) {
        decision_f = decision_function(test_chs[i], bias, lambdas, tags, learn_chs, learn_dims, kernel_char);
        if (decision_f * correct_sign > 0) correct++;
        fprintf(f, "%.15e\n", decision_f);
        printf("progress: %5.2f%%\r", (double)(i+1)/(double)n_test_chs * 100.);
        fflush(stdout);
    }
    printf("tests complete and written!\n");
    return (Int2){ .x = correct, .y = n_test_chs - correct };
}

ConfusionMatrix evaluate_set(Chain* eval_chains, int* eval_tags, int n_eval, Chain* learn_chains, int* learn_tags, int n_learn, double* lambdas, char kernel_char, int target_class) {
    printf("starting evaluation...\n");
    ConfusionMatrix cm = {0, 0, 0, 0};
    double bias = compute_bias(lambdas, learn_tags, learn_chains, n_learn, kernel_char);
    
    char time_str[80];
    time_t now = time(NULL);
    strftime(time_str, sizeof(time_str), "%Y-%m-%d_%H:%M:%S", localtime(&now));
    char run_dir[MAX_STR_LEN];
    sprintf(run_dir, "results/svm/decision-function-tests/kernel-%c/%s", kernel_char, time_str);
    mkdir_p(run_dir);
    printf("directory '%s' created!\n", time_str);

    char filename[2*MAX_STR_LEN];
    sprintf(filename, "%s/test-results.txt", run_dir);

    FILE *f = get_file(filename, "w");
    
    for (int i = 0; i < n_eval; i++) {
        double dec = decision_function(eval_chains[i], bias, lambdas, learn_tags, learn_chains, n_learn, kernel_char);
        int predicted = (dec > 0) ? 1 : -1;
        fprintf(f, "%d\t%.15e\t%.15e\n", eval_tags[i], dec, dec*eval_tags[i]);
        
        if (eval_tags[i] == target_class) {
            if (predicted == target_class) cm.TP++;
            else cm.FN++;
        } else {
            if (predicted == target_class) cm.FP++;
            else cm.TN++;
        }
        printf("progress: %5.2f%%\r", (double)(i+1)/(double)n_eval * 100.);
        fflush(stdout);
    }
    printf("set evaluation complete!\n\n");
    return cm;
}