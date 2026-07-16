#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <omp.h>

#ifdef _WIN32
    #include <direct.h>
    #define MKDIR(path) _mkdir(path)
#else
    #include <sys/stat.h>
    #include <sys/types.h>
    #define MKDIR(path) mkdir(path, 0755)
#endif

/*  lines in each seq file:
learn human:    1309
learn mouse:    373
test human:     1388
test mouse:     1379
*/

//**********PARALLEL THREADS**********
#define MAX_WORKER_PERCENTAGE 75.0
//************************************

//**********PARISI RAPUANO*************
#define NormRANu (2.3283063671E-10F)

extern unsigned int irr[256];
extern unsigned int ir1;
extern unsigned char ind_ran,ig1,ig2,ig3;

extern float Random(void);
extern void ini_ran(int SEMILLA);
//************************************



/*Global variables*/
#define N_AACIDS 21
#define N_PROPERTIES 9
#define CHAINLEN 298
#define EPSILON 1e-10 //Avoid division by 0
#define MAX_STR_LEN 256 //Max_str_len will really be thrice this value
/*Energy penalizations & weights*/
extern int PROPS_AA[N_AACIDS][N_PROPERTIES];
extern int CHAR_TO_INT_LUT[256];
#define ZERO_FREQ_PENALTY_LOG 100000 //Energy to sum for a zero-frequency AA in log humanness energy
#define ZERO_FREQ_PENALTY_LINEAR 6 //Energy to sum for a zero-frequency AA in linear humanness energy
#define ZERO_FREQ_PENALTY_PROPERTIES_DISTANCE 6700 //Energy to sum for a zero-frequency AA in properties distance energy
#define SVM_DISTANCE_PROP_FACTOR 1 // idk for the moment
#define SVM_PARAMETER_LIMIT 1

/*Metropolis parameters*/
extern double AA_MUTATION_PENALTY;
#define LAMBDA 3.0
#define STARTING_TARGET_ACCEPTANCE 0.25
#define MAX_AMORTIG 4.0
extern double WEIGHT_LOG;
extern double WEIGHT_PROP;
extern double WEIGHT_PENALTY;

/*Files*/
#define SEQS "seqs/"
#define RESULTS "results/"
#define METROPOLIS "metropolis/"
#define TXT ".txt"

#define FILE_L_MOUSE "learn_mouse"
#define L_MOUSE_N_LINES 373
#define FILE_L_HUMAN "learn_human"
#define L_HUMAN_N_LINES 1309
#define FRECS "_freqs"
#define ENTROPIS "_entropies"
#define CORREL "_correlations"

#define DECIMAL_PRECISION 5
#define BETA_PRECISION 10
#define ACCEPTANCE_PRECISION 3
#define ENERGY_PRECISION 8

/*Aminoacid and properties definition*/
//Order of properties: "-", HYDROPHOBIC, AROMATIC, ALIPHATIC, POLAR, SMALL, MINUSCULE, CHARGEDPLUS, CHARGEDMINUS -- note that the "-" is necessary, or else entropies would return -infty
#define HYDROPHOBIC "ACFGHIKLMTVWY"
#define AROMATIC "FHWY"
#define ALIPHATIC "ILV"
#define POLAR "DEHKNQRSTWY"
#define SMALL "ACDGNPSTV"
#define MINUSCULE "AGS"
#define CHARGEDPLUS "HKR"
#define CHARGEDMINUS "DE"

extern char AMINOACIDS[N_AACIDS + 1];
extern char *PROPERTIES[N_PROPERTIES];

/*Enums:*/


/*Structs:*/

// -  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
// HYDROPHOBIC
// AROMATIC
// ALIPHATIC
// POLAR
// SMALL
// MINUSCULE
// CHARGEDPLUS
// CHARGEDMINUS
typedef struct {
    double elements[N_AACIDS];
    double log_elements[N_AACIDS];//we can optimize calculations in metropolis_sweep calling log only once
    double properties[N_PROPERTIES];
} Aacid;

typedef struct {
    Aacid aas[CHAINLEN];
} Chain;

typedef struct {
    double x, y;
} Vec2;

typedef struct {
    double saa, laa, raa, taa;
    double spp, lpp, rpp, tpp;
    // Shannon AA, linear AA, Renyi AA, Tsallis AA, 
    // Shannon Props, linear Props, Renyi Props, Tsallis Props
} Entropies;

typedef struct {
    double log_humanness;
    double property_distance;
    double wanderer_penalty;
} Energy;



/*Functions:*/

//utils.c
void med_var(double* data, double* mean, double* variance, int n);
void minmax(double* data, double* max, double* min, int n);
int negative_chain(const Chain* ch);
FILE *get_file(char* filename, char* mode);
void sort_array(double* array, int n);

//metropolis.c
void generate_betas(double** betas, int n_betas, double* entropies, int n_entropies, double scale_factor, double epsilon, double cooling_rate);
void metropolis_sweep(char* murine_seq, const int* original_murine_indices, const Chain* human_ref_seq, double* local_beta, 
    double* acceptance, int chainlen, double w_log, double w_prop, double w_penalty);
int run_metropolis(char* murine_seq, const Chain* human_ref_seq, int n_sweeps, double** betas, int n_betas, char* filename);
void mega_metropolis(char* murine_seeds_filename, char* human_filename, int n_human_lines, int n_sweeps, double** betas, int n_betas, int n_metropolis);

//parsing.c
int char_to_int(char X);
char int_to_char(int X);
void initialize_char_to_int_LUT();
void initialize_properties_matrix();
void get_next_chain(FILE *f, Chain* out);
void append_file_to_chain_vector(char* filename, int n_lines, Chain *output, int starting_idx);
void read_next_line(FILE *f, char* out);

//chain-operations.c
void aacid_direct_sum(const Aacid* a, const Aacid* b, Aacid* out);
void aa_scale_only_aacids(const Aacid* aa, double scalar, Aacid* out);
void aa_scale_only_properties(const Aacid* aa, double scalar, Aacid* out);
void aa_normalize_properties(const Aacid* aa, Aacid* out);
void aa_normalize_aacids(const Aacid* aa, Aacid* out);
void chain_direct_sum(const Chain* a, const Chain* b, Chain* out);
void chain_deep_copy(const Chain* a, Chain* b);
void ch_normalize(const Chain* c, Chain* out);
void file_megaAacids(char *filename, int n_lines, Chain* out);
void entropy_vector(const Chain* mega_chain, Vec2 *output, char type, double order);
void all_entropies(const Chain* mega_chain, Entropies *output, double order);

void print_chain(const Chain* c);
void print_chain_to_file(const Chain* c, char* filename);
void print_entropies(Entropies *S);
void print_entropies_to_file(Entropies *S, char* filename);
/*WIP*/
void weigh_entropies(Entropies* input, double* output, double weighs[8]);

//correlations.c
void aa_correlation_matrix(Chain *input_chains, int n_chains, double output[CHAINLEN][CHAINLEN]);
void print_correlation_matrix_to_file(FILE *f, double corr[CHAINLEN][CHAINLEN]);
double scalar_product(double *a, double *b, int dims);
void correlation_for_position_and_aacid(Chain* chains, int n_chains, int position_index, int aa_index, Chain* out);
void correlation_for_position_and_prop(Chain* chains, int n_chains, int position_index, int prop_index, Chain* out);
void first_order_correlations(Chain* input_chains, int n_chains, Chain out_aacid[][N_AACIDS], Chain out_prop[][N_PROPERTIES]);