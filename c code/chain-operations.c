#include "head.h"

// sums frequencies of two aa
Aacid aacid_direct_sum(Aacid a, Aacid b) {
    Aacid out;
    for (int i = 0; i < N_AACIDS; i++) {
        out.elmts[i] = a.elmts[i] + b.elmts[i];
    }
    return out;
}

// scales frequencies of an aa
Aacid aa_freq_scale(Aacid aa, double scalar) {
    Aacid out;
    for (int i = 0; i < N_AACIDS; i++) {
        out.elmts[i] = aa.elmts[i] * scalar;
    }
    return out;
}

// sums frequencies of two chains
Chain chain_direct_sum(Chain a, Chain b) {
    Chain out;
    for (int i = 0; i < CHAINLEN; i++) {
        out.aas[i] = aacid_direct_sum(a.aas[i], b.aas[i]);
    }
    return out;
}

Chain ch_freq_scale(Chain c, double scalar) {
    Chain out;
    for (int i = 0; i < CHAINLEN; i++) {
        out.aas[i] = aa_freq_scale(c.aas[i], scalar);
    }
    return out;
}

// returns the schrödinger chain of a file, needs file and number of lines
Chain file_megaAacids(FILE *f, int n_lines) {
    Chain out = get_nex_chain(f);
    for (int i = 1; i < n_lines; i++) {
        out = chain_direct_sum(out, get_nex_chain(f));
    }
    return ch_freq_scale(out, 1/(double)n_lines);
}

double aa_shannon_entropy(Aacid aa) {
    for (int i = 0; i < N_AACIDS; i++) {
        
    }
}

void shannon_entropy_vector(Chain mega_chain, double *output) {

}