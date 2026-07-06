#include "head.h"

double scalar_product(double *a, double *b, int dims) {
    double sum = 0.;
    for (int i = 0; i < dims; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

CorrVec correlation_for_position(Chain chain, int position_index) {

}