#include "head.h"

double scalar_product(double *a, double *b, int dims) {
    double sum = 0.;
    for (int i = 0; i < dims; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

// takes many chains of 1's and 0's (not frequencies), a position and an aa to check the conditioned frequencies given that aa in that position
// the aminoacid to check for conditioning must be entered as the index in the 
Chain correlation_for_position_and_specific_aacid(Chain* chains, int n_chains, int position_index, int aa_index) {
    Chain out = chains[0];
    for (int ch_idx = 1; ch_idx < n_chains; ch_idx++) {
        if (chains[ch_idx].aas[position_index].elmts[aa_index] == 1.) {
            out = chain_direct_sum(out, chains[ch_idx]);
            if (ch_idx%100 == 0) printf("%d read lines\n", ch_idx);
        }
    }
    return ch_normalize(out);
}

// takes many chains of 1's and 0's (not frequencies), a position and a property to check the conditioned frequencies given that property in that position
Chain correlation_for_position_and_specific_aacid(Chain* chains, int n_chains, int position_index, int aa_index) {
    Chain out = chains[0];
    for (int ch_idx = 1; ch_idx < n_chains; ch_idx++) {
        if (chains[ch_idx].aas[position_index].elmts[aa_index] == 1.) {
            out = chain_direct_sum(out, chains[ch_idx]);
            if (ch_idx%100 == 0) printf("%d read lines\n", ch_idx);
        }
    }
    return ch_normalize(out);
}
