#include "head.h"

void med_var(double* data, double* mean, double* variance, int n) {
    int i;
    double sum, res;
    *mean=0;
    *variance=0;

    sum = 0.0;
    for (i=0; i<n; i++){
        sum += data[i];
    }
    *mean = sum/n;

    res = 0.0;
    for(i=0; i<n; i++){
        res += (data[i] - *mean) * (data[i] - *mean);//pow makes e^ln
    }

    *variance = res/(n-1);
}

void minmax(double* data, double* max, double* min, int n) {
    *max = *min = *data;
    for(int i=0; i<n; i++) {
        if(*max < *(data+i)) *max = *(data+i);
        if(*min > *(data+i)) *min = *(data+i);
    }
}

FILE *get_file(char* filename, char* mode) {
    FILE* f = fopen(filename, mode);
    if (f == NULL) fprintf(stderr, "Error: Could not open %s\n", filename);
    return f;
}