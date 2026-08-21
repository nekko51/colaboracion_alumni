#include "head.h"
/**
 * 's' for sigmoid (tanh)
 * 'p' for polynomial
 * 'l' for linear
 * 't' for tanimoto
 * 'c' for cauchy
 * 'r' for rational quadratic
 * 'g' for gaussian (RBF)
 * 'w' for piecewise function
 */
#define ACTUAL_KERNEL 'w'
#define PW_CENTER .2
#define PW_WIDTH .2

double AA_MUTATION_PENALTY;
double WEIGHT_LOG;
double WEIGHT_PROP;
double WEIGHT_PENALTY;

void initialize();
void svm_initNgram(char kernel_char);
void svm_initNannealing(double epsilon, int its_each_beta, int n_betas, int initial_beta_its, int print_svm_to_file, char kernel_char, double C_limit);
void svm_many_tests(double *C_candidates, int num_C, int runs_per_C, int n_betas, int its_for_initial_beta, int its_per_beta_annealing, char actual_kernel, double pw_center, double pw_width);

int main() {
    double C_candidates[] = {0.1, 1.0, 10.0, 100.0};
    int num_C = 4;
    int runs_per_C = 10;
    int its_for_initial_beta = 10000;
    int its_per_beta_annealing = 10000;
    int n_betas = 50;
    double center = .9;
    double width = .2;


        svm_many_tests(C_candidates, num_C, runs_per_C, n_betas, its_for_initial_beta, its_per_beta_annealing, 'w', center, width);



}




void initialize() {
    ini_ran(time(NULL));
    initialize_properties_matrix();
    initialize_char_to_int_LUT();
}

void svm_initNgram(char kernel_char) {
    printf("\nstarting initialization...\n");

    initialize();

    /*gram matrix allocation*/
    double **gram = malloc(N_DATA_POINTS * sizeof(double*));
    for (int i = 0; i < N_DATA_POINTS; i++) {
        gram[i] = malloc(N_DATA_POINTS * sizeof(double));
    }

    /*chains for gram matrix*/
    Chain *chs = malloc(N_DATA_POINTS * sizeof(Chain));
    append_file_to_chain_vector(SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES, chs, 0);
    append_file_to_chain_vector(SEQS FILE_L_MOUSE TXT, L_MOUSE_N_LINES, chs, L_HUMAN_N_LINES);
    printf("initialization complete!\n\n");

    svm_gram_matrix(chs, N_DATA_POINTS, gram, kernel_char, PW_CENTER, PW_WIDTH);
    char result_gram[MAX_STR_LEN];
    switch (kernel_char) {
    case 's':   
        sprintf(result_gram, "results/svm/gram/gram_matrix_sigmoid.txt");
        break;
    case 'p':   
        sprintf(result_gram, "results/svm/gram/gram_matrix_poly_o%d.txt", SVM_KERNEL_POW);
        break;
    default:
        sprintf(result_gram, "results/svm/gram/gram_matrix_gaussian.txt");
        break;
    };

    print_matrix_to_file(gram, N_DATA_POINTS, result_gram);

    for (int i = 0; i < N_DATA_POINTS; i++) free(gram[i]);
    free(gram);
    free(chs);
}



/**
 * - epsilon:           10⁻² ~ 10⁰
 * - initial_beta_its:  10²  ~ 10⁴
 */
void svm_initNannealing(double epsilon, int its_each_beta, int n_betas, int initial_beta_its, int print_svm_to_file, char kernel_char, double C_limit) {
    printf("\nstarting initialization...\n");

    ini_ran(time(NULL));
    initialize_properties_matrix();
    initialize_char_to_int_LUT();

    /*gram matrix allocation*/
    double **gram = malloc(N_DATA_POINTS * sizeof(double*));
    for (int i = 0; i < N_DATA_POINTS; i++) {
        gram[i] = malloc(N_DATA_POINTS * sizeof(double));
    }

    /*chain signs*/
    int *chs_s = malloc(N_DATA_POINTS * sizeof(int));
    for (int i = 0; i < L_HUMAN_N_LINES; i++) chs_s[i] = 1;
    for (int i = 0; i < L_MOUSE_N_LINES; i++) chs_s[i + L_HUMAN_N_LINES] = -1;


    char gram_file[MAX_STR_LEN];
    switch (kernel_char) {
    case 's':   
        sprintf(gram_file, "results/svm/gram/gram_matrix_sigmoid.txt");
        break;
    case 'p':   
        sprintf(gram_file, "results/svm/gram/gram_matrix_poly_o%d.txt", SVM_KERNEL_POW);
        break;
    default:
        sprintf(gram_file, "results/svm/gram/gram_matrix_gaussian.txt");
        break;
    };

    get_matrix_from_file(gram_file, N_DATA_POINTS, gram);

    double *lambdas = calloc(N_DATA_POINTS, sizeof(double));
    svm_initialize_lambdas(lambdas, chs_s, N_DATA_POINTS, C_limit);

    double *bets = malloc(n_betas * sizeof(double));
    bets[0] = svm_initial_beta(lambdas, epsilon, initial_beta_its, gram, chs_s, N_DATA_POINTS, C_limit);
    double cooling_rate = 1.1; 
    for (int i = 1; i < n_betas; i++) bets[i] = bets[i-1] * cooling_rate;

    printf("initialization complete!\n\n");

    svm_annealing(gram, chs_s, N_DATA_POINTS, epsilon, lambdas, bets, n_betas, its_each_beta, print_svm_to_file, C_limit, 0, kernel_char);

    free(chs_s);
    free(lambdas);
    for (int i = 0; i < N_DATA_POINTS; i++) free(gram[i]);
    free(gram);
}


// double C_candidates[] = {0.1, 1.0, 10.0, 100.0};
// int num_C = 4;
// int runs_per_C = 10;
// int its_for_initial_beta = 10000;
// int its_per_beta_annealing = 10000;
void svm_many_tests(double *C_candidates, int num_C, int runs_per_C, int n_betas, int its_for_initial_beta, int its_per_beta_annealing, char actual_kernel, double pw_center, double pw_width) {
    initialize();

    Chain *learn_chs, *val_chs, *test_chs;
    int *learn_tags, *val_tags, *test_tags;
    
    int n_learn = load_labeled_dataset("seqs/learn.txt", &learn_chs, &learn_tags);
    int n_val = load_labeled_dataset("seqs/validation.txt", &val_chs, &val_tags);
    int n_test = load_labeled_dataset("seqs/test.txt", &test_chs, &test_tags);

    if (n_learn == 0 || n_val == 0 || n_test == 0) {
        fprintf(stderr, "failed to read datasets\n");
    }

    double **gram_learn = malloc(n_learn * sizeof(double*));
    for (int i = 0; i < n_learn; i++) gram_learn[i] = malloc(n_learn * sizeof(double));

    char matrix_filename[MAX_STR_LEN];
    if (actual_kernel == 'w') {
        sprintf(matrix_filename, "results/svm/gram/piecewise/gram_%c_C-%.2lf_W-%.2lf.txt", actual_kernel, pw_center, pw_width);
    } else {
        sprintf(matrix_filename, "results/svm/gram/gram_%c.txt", actual_kernel);
    }
    FILE *file_check = fopen(matrix_filename, "r");
    if (file_check) {
        fclose(file_check);
        get_matrix_from_file(matrix_filename, n_learn, gram_learn);
    } else {
        svm_gram_matrix(learn_chs, n_learn, gram_learn, actual_kernel, pw_center, pw_width);
        print_matrix_to_file(gram_learn, n_learn, matrix_filename);
    }
    
    double best_val_f1 = -1.0;
    double best_C = 0.0;
    double *best_overall_lambdas = calloc(n_learn, sizeof(double));

    for (int c_idx = 0; c_idx < num_C; c_idx++) {
        double current_C = C_candidates[c_idx];
        double best_energy = 1e30;
        double *current_best_lambdas = calloc(n_learn, sizeof(double));

        #pragma omp parallel for
        for (int run = 0; run < runs_per_C; run++) {
            double *lambdas = calloc(n_learn, sizeof(double));
            svm_initialize_lambdas(lambdas, learn_tags, n_learn, current_C);
            
            double epsilon = 0.1; 
            int n_betas = 50;
            double *bets = malloc(n_betas * sizeof(double));
            bets[0] = svm_initial_beta(lambdas, epsilon, its_for_initial_beta, gram_learn, learn_tags, n_learn, current_C);
            for (int i = 1; i < n_betas; i++) bets[i] = bets[i-1] * 1.1;

            svm_annealing(gram_learn, learn_tags, n_learn, epsilon, lambdas, bets, n_betas, its_per_beta_annealing, 1, current_C, run, actual_kernel);
            
            double e = svm_energy(gram_learn, learn_tags, lambdas, n_learn);
            
            #pragma omp critical
            {
                if (e < best_energy) {
                    best_energy = e;
                    memcpy(current_best_lambdas, lambdas, n_learn * sizeof(double));
                }
            }
            
            free(lambdas);
            free(bets);
        }

        ConfusionMatrix val_cm = evaluate_set(val_chs, val_tags, n_val, learn_chs, learn_tags, n_learn, current_best_lambdas, actual_kernel, 1, pw_center, pw_width);
        double f1_val = 2.0 * val_cm.TP / (2.0 * val_cm.TP + val_cm.FP + val_cm.FN);

        if (f1_val > best_val_f1) {
            best_val_f1 = f1_val;
            best_C = current_C;
            memcpy(best_overall_lambdas, current_best_lambdas, n_learn * sizeof(double));
        }
        
        free(current_best_lambdas);
    }

    ConfusionMatrix test_cm_h = evaluate_set(test_chs, test_tags, n_test, learn_chs, learn_tags, n_learn, best_overall_lambdas, actual_kernel, 1, pw_center, pw_width);
    ConfusionMatrix test_cm_m = evaluate_set(test_chs, test_tags, n_test, learn_chs, learn_tags, n_learn, best_overall_lambdas, actual_kernel, -1, pw_center, pw_width);

    double F1H = 2.0 * test_cm_h.TP / (2.0 * test_cm_h.TP + test_cm_h.FP + test_cm_h.FN);
    double F1M = 2.0 * test_cm_m.TP / (2.0 * test_cm_m.TP + test_cm_m.FP + test_cm_m.FN);

    printf("\noptimal results (C = %.2f, kernel:%c)\n", best_C, actual_kernel);
    printf("H: TP:%4d, FP:%4d, TN:%4d, FN:%4d | F1H: %.15f\n", test_cm_h.TP, test_cm_h.FP, test_cm_h.TN, test_cm_h.FN, F1H);
    printf("M: TP:%4d, FP:%4d, TN:%4d, FN:%4d | F1M: %.15f\n", test_cm_m.TP, test_cm_m.FP, test_cm_m.TN, test_cm_m.FN, F1M);

    free(best_overall_lambdas);
    for (int i = 0; i < n_learn; i++) free(gram_learn[i]);
    free(gram_learn);
    free(learn_chs); free(learn_tags);
    free(val_chs); free(val_tags);
    free(test_chs); free(test_tags);
}