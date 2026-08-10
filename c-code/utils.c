#include "head.h"

/**
 * @brief Recursively creates a directory path.
 *
 * @param path The full directory path to create.
 * @return 0 on success, -1 on error.
 */
int mkdir_p(const char *path) {
    char *p;
    char *temp = strdup(path);
    if (temp == NULL) {
        fprintf(stderr, "Error: strdup failed in mkdir_p\n");
        return -1;
    }

    for (p = temp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            MKDIR(temp);
            *p = '/';
        }
    }
    MKDIR(temp);
    free(temp);
    return 0;
}

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

//returns 1 (true) if the chain is invalid, 0 (false) if it is invalid
int negative_chain(const Chain* ch) {
    if(ch->aas->elements[0] < -0.5) {
        fprintf(stderr, "Error: Chain returned negative values (%lf); stopping execution...\n", ch->aas->elements[0]);
        return 1;
    }
    else return 0;
}


FILE *get_file(char* filename, char* mode) {
    FILE* f = fopen(filename, mode);
    if (f == NULL) fprintf(stderr, "Error: Could not open %s\n", filename);
    return f;
}

int isSorted(double* array, int n) {
    for(int i=0; i<n; i++) {
        if(array[i]>array[i+1]) {
            return(0);
        }
    }
    return(1);
}

double* cosmicraysort(double* array, int n) {
    while(1) {//wait for cosmic rays to flip bits and sort our array
        if(isSorted(array, n)) {
            return(array);
        }
    }
}

int comp(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

void sort_array(double* array, int n) {
    qsort(array, n, sizeof(double), comp);
}