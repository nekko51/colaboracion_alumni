#include "head.h"

/*  lines in each seq file:
learn human:    1309
learn mouse:    373
test human:     1388
test mouse:     1379
*/

/** 
 * @brief file definitions in head.h
 * 
 * @param SEQS "seqs/"
 * @param RESULTS "results/"
 * @param TXT
 * @param FILE_L_MOUSE "learn_mouse"
 * @param L_MOUSE_N_LINES 373
 * @param FILE_L_HUMAN "learn_human"
 * @param L_HUMAN_N_LINES 1309
 * @param FRECS "_freqs"
 * @param ENTROPIS "_entropies"
 *
*/

double AA_MUTATION_PENALTY;
double WEIGHT_LOG;
double WEIGHT_PROP;
double WEIGHT_PENALTY;

void initialize(Chain* human_ref) {
    ini_ran(time(NULL));
    initialize_properties_matrix();
    initialize_char_to_int_LUT();
    file_megaAacids(SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES, human_ref);
}

// int main() {//i'm so happy i don't have to free every single malloc'd array if there's an error; if we had to, I think it's the only use of GOTO that wouldn't get you fired
//     /*Initial parameters*/
//     AA_MUTATION_PENALTY = 10;//0.5 would be less than almost every human mutation; 2.0 if it's clearly more human; 4.0 is pretty conservative (used log values for this)
//     //weight of 0.35, 0.35, 0.3 & penalty of 4 => 57.99 avg of avg hamming dist; penalty of 10 => 18.02 avg of avgs
//     WEIGHT_LOG = 0.35;//weights must sum to 1
//     WEIGHT_PROP = 0.35;
//     WEIGHT_PENALTY = 0.3;
//     int n_betas = 50;
//     int n_sweeps = 800;
//     int n_metropolis = 20;
//     int n_entropies = CHAINLEN;
//     double entropy_order_q = 0.34;
//     double scale_factor = 1.0;
//     double cooling_rate = 1.02;
//     double weighs[8] = {0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0};
//     Chain human_ref;

//     /*Variable & RNG Initializations*/
//     initialize(&human_ref);
    
//     double** betas = malloc(n_betas*sizeof(double*));
//     if(betas == NULL) {
//         fprintf(stderr, "Couldn't assign memory to betas matrix\n");
//         return(1);
//     }
//     for(int i=0; i<n_betas; i++) {
//         betas[i] = malloc(CHAINLEN*sizeof(double));
//         if(betas[i] == NULL) {
//             fprintf(stderr, "Couldn't assign memory to betas column %d/%d\n", i+1, n_betas);
//             return(1);
//         }
//     }

//     Entropies* oll_entropies = malloc(CHAINLEN*sizeof(Entropies));
//     if(oll_entropies == NULL) {
//         fprintf(stderr, "Couldn't assign memory to all_entropies array\n");
//         return(1);
//     }
//     double* entropies = malloc(CHAINLEN*sizeof(double));
//     if(entropies == NULL) {
//         fprintf(stderr, "Couldn't assign memory to entropies array\n");
//         return(1);
//     }

//     all_entropies(&human_ref, oll_entropies, entropy_order_q);
//     weigh_entropies(oll_entropies, entropies, weighs);//unfinished function, currently returns 1/2*(saa+spp)
//     print_entropies(oll_entropies);
//     generate_betas(betas, n_betas, entropies, n_entropies, scale_factor, EPSILON, cooling_rate);
//     mega_metropolis(SEQS FILE_L_MOUSE TXT, SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES, n_sweeps, betas, n_betas, n_metropolis);




//     /*Code preparation*/
//     // Chain ch = file_megaAacids(SEQS FILE_L_MOUSE TXT, L_MOUSE_N_LINES);
//     // if(negative_chain(&ch) == 1) return 1;

//     // print_chain_to_file(ch, RESULTS FILE_L_MOUSE FRECS TXT);
    
//     // Entropies *entropy = malloc(CHAINLEN * sizeof(Entropies));
//     // all_entropies(ch, entropy, .5);

//     // print_entropies_to_file(entropy, RESULTS FILE_L_MOUSE ENTROPIS TXT);

//     /*Free memory*/
//     for(int i=0; i<n_betas; i++) {
//         free(betas[i]);
//     }
//     free(betas);
//     free(oll_entropies);
//     free(entropies);
    
//     return 0;
// }


int main() {

    printf("\nstarting initialization...\n");

    ini_ran(time(NULL));
    initialize_properties_matrix();
    initialize_char_to_int_LUT();

    double epsilon = 5e-4;
    int N_DATA_POINTS = L_HUMAN_N_LINES + L_MOUSE_N_LINES;

    /*chains for gram matrix*/
    // Chain *chs = malloc(N_DATA_POINTS * sizeof(Chain));
    // append_file_to_chain_vector(SEQS FILE_L_HUMAN TXT, L_HUMAN_N_LINES, chs, 0);
    // append_file_to_chain_vector(SEQS FILE_L_MOUSE TXT, L_MOUSE_N_LINES, chs, L_HUMAN_N_LINES);

    // svm_gram_matrix(chs, N_DATA_POINTS, gram);
    // print_matrix_to_file(gram, N_DATA_POINTS, "results/svm/gram/gram_matrix_gaussian.txt");

    /*chain signs*/
    int *chs_s = malloc(N_DATA_POINTS * sizeof(int));
    for (int i = 0; i < L_HUMAN_N_LINES; i++) chs_s[i] = 1;
    for (int i = 0; i < L_MOUSE_N_LINES; i++) chs_s[i + L_HUMAN_N_LINES] = -1;

    /*gram matrix allocation*/
    double **gram = malloc(N_DATA_POINTS * sizeof(double*));
    for (int i = 0; i < N_DATA_POINTS; i++) {
        gram[i] = malloc(N_DATA_POINTS * sizeof(double));
    }

    printf("starting gram matrix reading...\n");
    get_matrix_from_file("results/svm/gram/gram_matrix_gaussian.txt", N_DATA_POINTS, gram);
    printf("gram matrix reading complete!\n");

    double *lambdas = calloc(N_DATA_POINTS, sizeof(double));

    /*--------------------------------------------------------------------*/
    int its_each_beta = 100;
    int n_betas = 1000;
    /*--------------------------------------------------------------------*/

    double *bets = malloc(n_betas * sizeof(double));
    bets[0] = svm_initial_beta(lambdas, epsilon, 100, gram, chs_s, N_DATA_POINTS);
    double cooling_rate = 1.1; 
    for (int i = 1; i < n_betas; i++) bets[i] = bets[i-1] * cooling_rate;

    char time_str[80];
    time_t now = time(NULL);
    strftime(time_str, sizeof(time_str), "%Y-%m-%d_%H_%M_%S", localtime(&now));

    char run_dir[MAX_STR_LEN];
    sprintf(run_dir, "results/svm/annealing/%s_b%05d_its%05d", time_str, n_betas, its_each_beta);
    mkdir_p(run_dir);

    printf("initialization complete!\n\n");

    svm_annealing(gram, chs_s, N_DATA_POINTS, epsilon, lambdas, bets, n_betas, its_each_beta, 1, run_dir);

    // for(int i = 0; i < N_DATA_POINTS; i++) printf("%lf\n", lambdas[i]);


    free(chs_s);
    free(lambdas);
    for (int i = 0; i < N_DATA_POINTS; i++) free(gram[i]);
    free(gram);
    
    return 0;
}

/*  lines in each seq file:
learn human:    1309
learn mouse:    373
test human:     1388
test mouse:     1379
*/