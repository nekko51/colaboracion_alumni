#include "head.h"

void initialize();
void svm_initNgram(char kernel_char);
void svm_initNannealing(double epsilon, int its_each_beta, int n_betas, int initial_beta_its, int print_svm_to_file, char kernel_char);


int main() {
    svm_initNgram('p');
}

// int main() {

//     /*initialization*/
//     printf("\nstarting initialization...\n");

//     initialize();

//     //lambdas
//     double* lambdas = calloc(N_DATA_POINTS, sizeof(double));
//     get_best_lambdas_from_csv(lambdas, N_DATA_POINTS);

//     //tags
//     int *signs = malloc(N_DATA_POINTS * sizeof(int));
//     for (int i = 0; i < L_HUMAN_N_LINES; i++) signs[i] = 1;
//     for (int i = 0; i < L_MOUSE_N_LINES; i++) signs[i + L_HUMAN_N_LINES] = -1;

//     //gram matrix
//     double **gram = malloc(N_DATA_POINTS * sizeof(double*));
//     for (int i = 0; i < N_DATA_POINTS; i++) { gram[i] = malloc(N_DATA_POINTS * sizeof(double)); }
//     get_matrix_from_file("results/svm/gram/gram_matrix_sigmoid.txt", N_DATA_POINTS, gram);

//     //learn chains
//     Chain *learn_chs = malloc(N_DATA_POINTS * sizeof(Chain));
//     append_file_to_chain_vector(SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES, learn_chs, 0);
//     append_file_to_chain_vector(SEQS FILE_L_MOUSE TXT, L_MOUSE_N_LINES, learn_chs, L_HUMAN_N_LINES);

//     //test chains
//     Chain *test_h = malloc(T_HUMAN_N_LINES * sizeof(Chain));
//     append_file_to_chain_vector(SEQS FILE_T_HUMAN TXT, T_HUMAN_N_LINES, test_h, 0);
//     Chain *test_m = malloc(T_MOUSE_N_LINES * sizeof(Chain));
//     append_file_to_chain_vector(SEQS FILE_T_MOUSE TXT, T_MOUSE_N_LINES, test_m, 0);

//     printf("initialization complete!\n\n");
//     /**/

//     test_results_to_file(N_DATA_POINTS, lambdas, signs, learn_chs, test_h, T_HUMAN_N_LINES, 's');
//     test_results_to_file(N_DATA_POINTS, lambdas, signs, learn_chs, test_m, T_MOUSE_N_LINES, 's');

//     free(lambdas);
//     free(signs);
//     for (int i = 0; i < N_DATA_POINTS; i++) free(gram[i]);
//     free(gram);
//     free(learn_chs);
    
//     return 0;
// }

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

    svm_gram_matrix(chs, N_DATA_POINTS, gram, kernel_char);
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
void svm_initNannealing(double epsilon, int its_each_beta, int n_betas, int initial_beta_its, int print_svm_to_file, char kernel_char) {
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
    svm_initialize_lambdas(lambdas, chs_s, N_DATA_POINTS);

    double *bets = malloc(n_betas * sizeof(double));
    bets[0] = svm_initial_beta(lambdas, epsilon, initial_beta_its, gram, chs_s, N_DATA_POINTS);
    double cooling_rate = 1.1; 
    for (int i = 1; i < n_betas; i++) bets[i] = bets[i-1] * cooling_rate;

    char time_str[80];
    time_t now = time(NULL);
    strftime(time_str, sizeof(time_str), "%Y-%m-%d_%H:%M:%S", localtime(&now));
    char run_dir[MAX_STR_LEN];
    sprintf(run_dir, "results/svm/annealing/batches/betas%05d_iterations%05d/%s", n_betas, its_each_beta, time_str);
    mkdir_p(run_dir);
    printf("directory '%s' created!\n", time_str);

    printf("initialization complete!\n\n");

    svm_annealing(gram, chs_s, N_DATA_POINTS, epsilon, lambdas, bets, n_betas, its_each_beta, print_svm_to_file, run_dir);

    free(chs_s);
    free(lambdas);
    for (int i = 0; i < N_DATA_POINTS; i++) free(gram[i]);
    free(gram);
}