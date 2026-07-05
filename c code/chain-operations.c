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
Aacid aa_scale(Aacid aa, double scalar) {
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

// scales frequencies of a chain
Chain ch_scale(Chain c, double scalar) {
    Chain out;
    for (int i = 0; i < CHAINLEN; i++) {
        out.aas[i] = aa_scale(c.aas[i], scalar);
    }
    return out;
}

// returns the schrödinger chain of a file, needs file and number of lines
Chain file_megaAacids(char *filename, int n_lines) {
    FILE *f = fopen(filename, "r");
    if (f == NULL) {printf("no se pudo abrir %s\n", filename);return;}
    Chain out = get_nex_chain(f);
    for (int i = 1; i < n_lines; i++) {
        out = chain_direct_sum(out, get_nex_chain(f));
    }
    fclose(f);
    return ch_scale(out, 1/(double)n_lines);
}


/*entropies*/

double aa_shannon_entropy(Aacid aa) {
    double sum = 0;
    for (int i = 0; i < N_AACIDS; i++) {
        sum += -aa.elmts[i] * log(aa.elmts[i]);
    }
    return sum;
}

double aa_linear_entropy(Aacid aa) {
    double sum = 0;
    for (int i = 0; i < N_AACIDS; i++) {
        sum += -aa.elmts[i] * (1 - aa.elmts[i]);
    }
    return sum;
}

double aa_renyi_entropy(Aacid aa, double q) {
    double sum = 0;
    for (int i = 0; i < N_AACIDS; i++) {
        sum += pow(aa.elmts[i],q);
    }
    return log(sum)/(1-q);
}

double aa_tsallis_entropy(Aacid aa, double q) {
    double sum = 0;
    for (int i = 0; i < N_AACIDS; i++) {
        sum += pow(aa.elmts[i],q);
    }
    return (1 - sum) / ( q - 1);
}

/* calcs entropy for each entry in a chain, output must be at least N_AACIDS long
for argument ''type'':
- 'l':    linear entropy
- 'r':    Renyi entropy*
- 't':    Tsallis entropy*
- none of the above:  Shannon entropy (default)

[*] order \in [0,1), only used for these [*], otherwise ignored
*/
void entropy_vector(Chain mega_chain, double *output, char type, double order) {
    if (type == 'l') {
        //linear
        for (int i = 0; i < CHAINLEN; i++) {
            output[i] = aa_linear_entropy(mega_chain.aas[i]);
        }
    } else if (type == 'r') {
        //renyi
        for (int i = 0; i < CHAINLEN; i++) {
            output[i] = aa_renyi_entropy(mega_chain.aas[i], order);
        }
    } else if (type == 't') {
        //tsallis
        for (int i = 0; i < CHAINLEN; i++) {
            output[i] = aa_tsallis_entropy(mega_chain.aas[i], order);
        }
    } else {
        //shannon
        for (int i = 0; i < CHAINLEN; i++) {
            output[i] = aa_shannon_entropy(mega_chain.aas[i]);
        }
    }
}



