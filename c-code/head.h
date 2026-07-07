#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
/*  lines in each seq file:
learn human:    1309
learn mouse:    373
test human:     1388
test mouse:     1379
*/

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
#define MAX_STR_LEN 600
/*Energy penalizations & weights*/
extern int PROPS_AA[N_AACIDS][N_PROPERTIES];
#define ZERO_FREQ_PENALTY_LOG 100000 //Energy to sum for a zero-frequency AA in log humanness energy
#define ZERO_FREQ_PENALTY_LINEAR 6 //Energy to sum for a zero-frequency AA in linear humanness energy
#define ZERO_FREQ_PENALTY_PROPERTIES_DISTANCE 6700 //Energy to sum for a zero-frequency AA in properties distance energy

//must sum to 1
#define WEIGHT_LOG 0.5 
#define WEIGHT_PROP 0.5

/*Files*/
#define SEQS "seqs/"
#define RESULTS "results/"
#define TXT ".txt"

#define FILE_L_MOUSE "learn_mouse"
#define L_MOUSE_N_LINES 373
#define FILE_L_HUMAN "learn_human"
#define L_HUMAN_N_LINES 1309
#define FRECS "_freqs"
#define ENTROPIS "_entropies"

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
} Energy;



/*Functions:*/

//utils.c
void med_var(double* data, double* mean, double* variance, int n);
void minmax(double* data, double* max, double* min, int n);
int negative_chain(const Chain ch);
FILE *get_file(char* filename, char* mode);

//metropolis.c
double log_humanness_energy(const Chain* human_ref_seq, const char* seq, int n);
double linear_humanness_energy(const Chain* human_ref_seq, const char* seq, int n);
double property_distance_energy(const Chain* human_ref_seq, const char* seq, int n);
Energy energy_calculation(const Chain* ref, const char* seq, int n);
void metropolis_sweep(char* murine_seq, const Chain* human_ref_seq, double beta, double* acceptance, int n, double w_log, double w_prop);
int run_metropolis(char* murine_seq, int n_steps, double* betas, int n_betas);

//parsing.c
int char_to_int(char X);
char int_to_char(int X);
Chain get_next_chain(FILE *f);
void append_file_to_chain_vector(char* filename, int n_lines, Chain *output, int starting_idx);
void initialize_properties_matrix();

//chain-operations.c
Aacid aacid_direct_sum(Aacid a, Aacid b);
Aacid aa_scale_only_aacids(Aacid aa, double scalar);
Aacid aa_scale_only_properties(Aacid aa, double scalar);
Chain chain_direct_sum(Chain a, Chain b);
void chain_deep_copy(const Chain a, Chain* b);
Chain ch_normalize(Chain c);
Chain file_megaAacids(char *filename, int n_lines);
void entropy_vector(Chain mega_chain, Vec2 *output, char type, double order);
void all_entropies(Chain mega_chain, Entropies *output, double order);
void print_chain(Chain c);
void print_chain_to_file(Chain c, char* filename);
void print_entropies(Entropies *S);
void print_entropies_to_file(Entropies *S, char* filename);

//correlations.c
double scalar_product(double *a, double *b, int dims);
Chain correlation_for_position_and_aacid(Chain* chains, int n_chains, int position_index, int aa_index);
Chain correlation_for_position_and_prop(Chain* chains, int n_chains, int position_index, int prop_index);
void first_order_correlations(Chain* input_chains, int n_chains, Chain out_aacid[][N_AACIDS], Chain out_prop[][N_PROPERTIES]);