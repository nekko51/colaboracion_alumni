#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define fran 3

#define N_AACIDS 21
#define CHAINLEN 298

/*  lines in each seq file:
learn human:    1309
learn mouse:    373
test human:     1388
test mouse:     1379
*/

#define AACIDS "-ACDEFGHIKLMNPQRSTVWY"
#define HYDROPHOBIC "ACFGHIKLMTVWY"
#define AROMATIC "FHWY"
#define ALIPHATIC "ILV"
#define POLAR "DEHKNQRSTWY"
#define SMALL "ACDGNPSTV"
#define MINUSCULE "AGS"
#define CHARGEDPLUS "HKR"
#define CHARGEDMINUS "DE"

/*Structs:*/

// -  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
// 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
typedef struct {
    double elmts[N_AACIDS];
} Aacid;

typedef struct {
    Aacid aas[CHAINLEN];
} Chain;

/*Functions:*/

//utils.c
void med_var(double* data, double* mean, double* variance, int n);
void minmax(double* data, double* max, double* min, int n);

//metropolis.c
double calculate_energy(const char* seq1, const char* seq2, int length);
void metropolis_sweep(char* murine_seq, const char* human_ref, double beta, int n);

//parsing.c
Chain get_nex_chain(FILE *f);

