#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define fran 3

#define N_AACIDS 21
#define N_PROPERTIES 9
#define CHAINLEN 298

/*  lines in each seq file:
learn human:    1309
learn mouse:    373
test human:     1388
test mouse:     1379
*/

#define AACIDS "-ACDEFGHIKLMNPQRSTVWY"

//Order of properties:
//{"-", HYDROPHOBIC, AROMATIC, ALIPHATIC, POLAR, SMALL, MINUSCULE, CHARGEDPLUS, CHARGEDMINUS};
//Note that this "-" is necessary, or else entropies would return -infty
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
    double elmts[N_AACIDS];
    double props[N_PROPERTIES];
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

/*Functions:*/

//utils.c
void med_var(double* data, double* mean, double* variance, int n);
void minmax(double* data, double* max, double* min, int n);
FILE *get_file(char* filename, char* mode);

//metropolis.c
double calculate_energy(const char* seq1, const char* seq2, int length);
void metropolis_sweep(char* murine_seq, const char* human_ref, double beta, int n);

//parsing.c
int char_to_int(char X);
Chain get_next_chain(FILE *f);
void read_file(char* filename, int n_lines, Chain *output);

//chain-operations.c
Aacid aacid_direct_sum(Aacid a, Aacid b);
Aacid aa_scale_only_aacids(Aacid aa, double scalar);
Aacid aa_scale_only_props(Aacid aa, double scalar);
Chain chain_direct_sum(Chain a, Chain b);
Chain ch_normalize(Chain c);
Chain file_megaAacids(char *filename, int n_lines);
void entropy_vector(Chain mega_chain, Vec2 *output, char type, double order);
void all_entropies(Chain mega_chain, Entropies *output, double order);
void print_chain(Chain c);
void print_chain_to_file(Chain c, char* filename);
void print_entropies(Entropies *S);
void print_entropies_to_file(Entropies *S, char* filename);

//correlations.c
