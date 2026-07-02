#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define fran 3

/*Functions:*/

//utils.c
void med_var(double* data, double* mean, double* variance, int n);
void minmax(double* data, double* max, double* min, int n);

//metropolis.c
double calculate_energy(const char* seq1, const char* seq2, int length);
void metropolis_sweep(char* murine_seq, const char* human_ref, double beta, int n);