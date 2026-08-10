#include "head.h"

double aa_sq_distance(Aacid aa, Aacid ab) {
    double sum = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sum += (aa.elements[i] - ab.elements[i]) * (aa.elements[i] - ab.elements[i]);
    }
    return sum;
}

double chain_sq_distance(Chain a, Chain b) {
    double sum = 0.;
    for (int i = 0; i < CHAINLEN; i++) {
        sum += aa_sq_distance(a.aas[i], b.aas[i]);
    }
    return sum;
}

// @return exp( - SVD_DISTANCE_PROP_FACTOR * chain_sq_distance(a, b))
double custom_distance(Chain a, Chain b) {
    return exp( - SVM_DISTANCE_PROP_FACTOR * chain_sq_distance(a, b));
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

// @brief computes (signed) gram matrix for a series of learn chains
// @returns a matrix of size n_chains * n_chains
void svm_gram_matrix(Chain *chains, int n_chains, double **out) {
    printf("starting gram matrix computation...\n");
    for (int i = 0; i < n_chains; i++) {
        for (int j = 0; j < n_chains; j++) {
            out[i][j] = custom_distance(chains[i], chains[j]);
        }
        if (i%10 == 0)  printf("gram matrix progress:\t%5.2f%%\n", (double)(i+1)/(double)n_chains * 100);
    }
    printf("finished gram matrix!\n");
}

/**
 * decision function
 * f(x) = b + \sum_i \lambda_i s_i dist(x_i, x)
 * 
 * bias
 * b = s_p - \sum_i \lambda_i s_i dist(x_i, x_p) 
 * where \lambda_p > 0
 * 
 */

double compute_bias(double *lambda, double *signs, Chain *chains, int n_chains) {
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
        bias -= lambda[i] * signs[i] * custom_distance(chains[i], chains[idx]);
    }
    return bias;
}

double decision_function(Chain x_input, double bias, double *lambda, int *signs, Chain *chains, int n_chains) {
    double out = bias;
    for (int i = 0; i < n_chains; i++) {
        out += lambda[i] * signs[i] * custom_distance(x_input, chains[i]);
    }
    return out;
}


/**
 * maximize 
 * \sum_i \lambda_i - \frac{1}{2} \sum_i \sum_j \lambda_i \lambda_j s_i s_j * dist(x_i, x_j)
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

void svm_change(double *v, double *v_, int v_dim, double eps, int *signs) {
    do {
        v_[v_dim-1] = 0;
        for (int i = 0; i < v_dim - 1; i++) {
            do {
                v_[i] = v[i] + Random() * 2 * eps - eps;
            } while (v_[i] <= 0. || v_[i] >= SVM_PARAMETER_LIMIT);
            v_[v_dim-1] -= v_[i] * signs[i];
        }
        v_[v_dim-1] *= signs[v_dim-1]; 
    } while (v_[v_dim-1] <= 0. || v_[v_dim-1] >= SVM_PARAMETER_LIMIT);
}

int svm_metropolis(double e_old, double e_new, double beta) {
    double ce = exp(-beta*(e_new - e_old));
    if (ce > Random()) {
        return 1;
    }
    return 0;
}

double svm_initial_beta(double *initial_lambdas, double eps, int iterations, double **gram_matrix, int *signs, int dim) {
    printf("starting initial beta computation...\n");
    double *lambdas_ = malloc(dim * sizeof(double));
    double energy_old = svm_energy(gram_matrix, signs, initial_lambdas, dim);
    double sum = 0.0;
    
    for (int i = 0; i < iterations; i++) {
        svm_change(initial_lambdas, lambdas_, dim, eps, signs);
        double delta = svm_energy(gram_matrix, signs, lambdas_, dim) - energy_old;
        sum += fabs(delta);
        printf("progress:\t%5.2f%%\r", (double)(i+1)/(double)iterations * 100.);
        fflush(stdout);
    }
    
    free(lambdas_);
    printf("finished initial beta computation!\n");
    return (double)iterations / sum;
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
void svm_annealing(double **gram_matrix, int *signs, int dimension, double epsilon, double *lambdas, double *betas, int n_betas, int iterations_per_beta, int print_to_file, char *filename) {
    printf("starting simulated annealing...\n");
    double energy_old, energy_new; double *lambdas_ = calloc(dimension, sizeof(double)); int acceptance_count;
    energy_old = svm_energy(gram_matrix, signs, lambdas, dimension);
    printf("progress:\t%5.2f%%\r", 0.);
    fflush(stdout);
    if (print_to_file) {

        /*file things*/
        char summry_file_name[MAX_STR_LEN];
        char log_file_name[MAX_STR_LEN];
        sprintf(summry_file_name, "%s/summary.txt", filename);
        sprintf(log_file_name, "%s/full-log.txt", filename);

        FILE *sum_file = get_file(summry_file_name, "w");
        FILE *log_file = get_file(log_file_name, "w");

        fprintf(sum_file, "B-idx\tEnergy\tAcceptance\t");
        for (int i = 0; i < dimension; i++) fprintf(sum_file, "lambda_%d\t", i);
        fprintf(sum_file, "0\t%lf\t0.\t", energy_old);
        for (int i = 0; i < dimension; i++) fprintf(sum_file, "%lf\t", lambdas[i]);

        fprintf(log_file, "B-idx\tIt-idx\tEnergy\tAcceptance~\n");
        /*end of file things*/

        for (int beta_idx = 0; beta_idx < n_betas; beta_idx++) {
            acceptance_count = 0;

            for (int iteration_idx = 0; iteration_idx < iterations_per_beta; iteration_idx++) {
                svm_change(lambdas, lambdas_, dimension, epsilon, signs);
                energy_new = svm_energy(gram_matrix, signs, lambdas_, dimension);

                if (svm_metropolis(energy_old, energy_new, betas[beta_idx])) {
                    energy_old = energy_new;
                    acceptance_count++;
                    for (int i = 0; i < dimension; i++) {
                        lambdas[i] = lambdas_[i];
                    }
                }

                fprintf(log_file, "%d\t%d\t%g\t%g\n", beta_idx, iteration_idx, energy_old, (double)acceptance_count/iterations_per_beta);
            }

            fprintf(sum_file, "%d\t%lf\t", beta_idx+1, energy_old);
            fprintf(sum_file, "%lf\t", (double)acceptance_count/(double)iterations_per_beta);
            for (int i = 0; i < dimension; i++) fprintf(sum_file, "%lf\t", lambdas[i]);

            printf("progress:\t%5.2f%%\r", (double)(beta_idx+1)/(double)n_betas * 100.);
            fflush(stdout);

            if ((double)acceptance_count/(double)iterations_per_beta > 0.5) epsilon *= 1.05; 
            else                                                            epsilon *= 0.95;

        }
        fclose(sum_file);
        fclose(log_file);
    } else {
        for (int beta_idx = 0; beta_idx < n_betas; beta_idx++) {

            for (int iteration_idx = 0; iteration_idx < iterations_per_beta; iteration_idx++) {
                svm_change(lambdas, lambdas_, dimension, epsilon, signs);
                energy_new = svm_energy(gram_matrix, signs, lambdas_, dimension);

                if (svm_metropolis(energy_old, energy_new, betas[beta_idx])) {
                    energy_old = energy_new;
                    acceptance_count++;
                    for (int i = 0; i < dimension; i++) {
                        lambdas[i] = lambdas_[i];
                    }
                }
            }

            printf("progress:\t%5.2f%%\r", (double)(beta_idx+1)/(double)n_betas * 100.);
            fflush(stdout);

            if ((double)acceptance_count/(double)iterations_per_beta > 0.5) epsilon *= 1.05; 
            else                                                            epsilon *= 0.95;
        }
    }
    free(lambdas_);
    printf("finished simulated annealing!\n");
}
