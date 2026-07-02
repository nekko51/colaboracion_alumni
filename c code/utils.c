#include "head.h"

void med_var(double* data, int n, double* mean, double* variance) {
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