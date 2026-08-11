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

double kernel(Chain a, Chain b) {
    return gaussian_radial_basis_function(a, b);
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
            out[i][j] = kernel(chains[i], chains[j]);

        }
        printf("progress: %5.2f%%\r", (double)(i+1)/(double)n_chains * 100.);
        fflush(stdout);
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
        bias -= lambda[i] * signs[i] * kernel(chains[i], chains[idx]);
    }
    return bias;
}

double decision_function(Chain x_input, double bias, double *lambda, int *signs, Chain *chains, int n_chains) {
    double out = bias;
    for (int i = 0; i < n_chains; i++) {
        out += lambda[i] * signs[i] * kernel(x_input, chains[i]);
    }
    return out;
}


SvmMove svm_propose_move(double *v, int v_dim, double eps, int *signs) {
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
        double max_di = SVM_PARAMETER_LIMIT - v[move.i];
        double min_dj, max_dj;

        if (signs[move.i] == signs[move.j]) {
            min_dj = v[move.j] - SVM_PARAMETER_LIMIT;
            max_dj = v[move.j];
        } else {
            min_dj = -v[move.j];
            max_dj = SVM_PARAMETER_LIMIT - v[move.j];
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
        double max_di = SVM_PARAMETER_LIMIT - v[i];
        double min_dj, max_dj;

        if (signs[i] == signs[j]) {
            min_dj = v[j] - SVM_PARAMETER_LIMIT;
            max_dj = v[j];
        } else {
            min_dj = -v[j];
            max_dj = SVM_PARAMETER_LIMIT - v[j];
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

double svm_initial_beta(double *initial_lambdas, double eps, int iterations, double **gram_matrix, int *signs, int dim) {
    printf("starting initial beta computation...\n");
    
    double *F = calloc(dim, sizeof(double));
    for (int k = 0; k < dim; k++) {
        for (int m = 0; m < dim; m++) {
            F[k] += initial_lambdas[m] * signs[m] * gram_matrix[k][m];
        }
    }

    double sum = 0.0;
    
    for (int it = 0; it < iterations; it++) {
        SvmMove move = svm_propose_move(initial_lambdas, dim, eps, signs);
        
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
    
    printf("finished initial beta computation! b = %.3g\n", ret);
    return ret;
}

void svm_initialize_lambdas(double *lambdas, int *signs, int dim) {
    double sum_pos = 0.0;
    double sum_neg = 0.0;

    for (int i = 0; i < dim; i++) {
        lambdas[i] = Random() * (SVM_PARAMETER_LIMIT * 0.1);
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
        char summry_file_name[MAX_STR_LEN];
        char log_file_name[MAX_STR_LEN];
        sprintf(summry_file_name, "%s/summary.txt", filename);
        sprintf(log_file_name, "%s/full-log.txt", filename);

        FILE *sum_file = get_file(summry_file_name, "w");
        FILE *log_file = get_file(log_file_name, "w");

        fprintf(sum_file, "B-idx\tBeta\tEnergy\tAcceptance\n");

        fprintf(log_file, "Line\tB-idx\tBeta\tIt-idx\tEnergy\tAcceptance~\n");
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

                printf("progress:\t%5.2f%%\r", (double)(beta_idx*iterations_per_beta+iteration_idx+1)/(double)(n_betas*iterations_per_beta) * 100.);
                fflush(stdout);

                fprintf(log_file, "%d\t%d\t%lf\t%d\t%lf\t%g\n", iteration_idx+beta_idx*iterations_per_beta, beta_idx, betas[beta_idx], iteration_idx, energy_old, (double)acceptance_count/iterations_per_beta);
            }

            fprintf(sum_file, "%d\t%lf\t%lf\t%lf\n", beta_idx+1, betas[beta_idx], energy_old, (double)acceptance_count/(double)iterations_per_beta);

            if ((double)acceptance_count/(double)iterations_per_beta > 0.5) epsilon *= 1.05; 
            else                                                            epsilon *= 0.95;

        }

        char res_file_name[MAX_STR_LEN];
        sprintf(res_file_name, "%s/result-lambdas.txt", filename);
        FILE *res_file = get_file(res_file_name, "w");

        fprintf(res_file, "%g\t", energy_old);
        for(int i = 0; i < dimension; i++)  fprintf(res_file, "%g\t", lambdas[i]);

        fclose(res_file);
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

                printf("progress:\t%5.2f%%\r", (double)(beta_idx*iterations_per_beta+iteration_idx+1)/(double)(n_betas*iterations_per_beta) * 100.);
                fflush(stdout);
            }

            if ((double)acceptance_count/(double)iterations_per_beta > 0.5) epsilon *= 1.05; 
            else                                                            epsilon *= 0.95;
        }
    }
    free(lambdas_);
    printf("finished simulated annealing!\n");
}
