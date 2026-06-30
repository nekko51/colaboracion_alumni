#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define fran 3

/*Functions:*/
double calculate_energy(const char* seq1, const char* seq2, int length);

void metropolis_sweep(char* murine_seq, const char* human_ref, int n_steps, double beta, int n);