#include "head.h"


// @brief calculates correlation between two positions of la chain in a list of chains using the mutual information formula
// @param input_chains array of chains used as input
// @param n_chains number of chains in inputted array
// @param idx_a one of the indices in the chain
// @param idx_b the other index for correlation
// @returns statistical correlation between those two positions 
double aa_correlation(Chain *input_chains, int n_chains, int idx_a, int idx_b) {
    int count_a[N_AACIDS]          = {0};
    int count_b[N_AACIDS]          = {0};
    int count_ab[N_AACIDS][N_AACIDS] = {0};
    int valid = 0;

    for (int chain = 0; chain < n_chains; chain++) {
        int val_a = -1, val_b = -1;
        for (int i = 0; i < N_AACIDS; i++) {
            if (input_chains[chain].aas[idx_a].elements[i] >= 1.0 - EPSILON) val_a = i;
            if (input_chains[chain].aas[idx_b].elements[i] >= 1.0 - EPSILON) val_b = i;
        }
        if (val_a != -1 && val_b != -1) {
            count_a[val_a]++;
            count_b[val_b]++;
            count_ab[val_a][val_b]++;
            valid++;
        }
    }

    if (valid == 0) return 0.0;

    // mutual information
    double mi = 0.0;
    for (int a = 0; a < N_AACIDS; a++) {
        for (int b = 0; b < N_AACIDS; b++) {
            if (count_ab[a][b] > 0) {
                double p_ab = (double)count_ab[a][b] / valid;
                double p_a  = (double)count_a[a]     / valid;
                double p_b  = (double)count_b[b]     / valid;
                mi += p_ab * log2(p_ab / (p_a * p_b));
            }
        }
    }

    // marginal entropies
    double h_a = 0.0, h_b = 0.0;
    for (int i = 0; i < N_AACIDS; i++) {
        if (count_a[i] > 0) {
            double p = (double)count_a[i] / valid;
            h_a -= p * log2(p);
        }
        if (count_b[i] > 0) {
            double p = (double)count_b[i] / valid;
            h_b -= p * log2(p);
        }
    }

    double denom = (h_a < h_b) ? h_a : h_b;
    if (denom < EPSILON) return 0.0;

    return mi / denom;
}

// @brief calculates correlation matrix
// @param input_chains array of chains used as input
// @param n_chains number of chains in inputted array
// @returns output matrix used as output
void aa_correlation_matrix(Chain *input_chains, int n_chains, double output[CHAINLEN][CHAINLEN]) {
    for (int a = 0; a < CHAINLEN; a++) {
        for (int b = 0; b < CHAINLEN; b++) {
            output[a][b] = aa_correlation(input_chains, n_chains, a, b);
        }
    }
}

void print_correlation_matrix_to_file(FILE *f, double corr[CHAINLEN][CHAINLEN]) {
    for (int i = 0; i < CHAINLEN; i++) {
        for (int j = 0; j < CHAINLEN; j++) {
            fprintf(f, "%lf\t", corr[i][j]);
        }
        fprintf(f, "\n");
    }
}






double scalar_product(double *a, double *b, int dims) {
    double sum = 0.;
    for (int i = 0; i < dims; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}



// ----------------- form here, not used code ---------------------


// takes many chains of 1's and 0's (not frequencies), a position and an aa to check the conditioned frequencies given that aa in that position
// the aminoacid to check for conditioning must be entered as the index in the AMINOACIDS list
void correlation_for_position_and_aacid(Chain* chains, int n_chains, int position_index, int aa_index, Chain* out) {
    memset(out, 0, sizeof(Chain));
    for (int ch_idx = 0; ch_idx < n_chains; ch_idx++) {
        if (chains[ch_idx].aas[position_index].elements[aa_index] >= 1.-EPSILON) {
            chain_direct_sum(out, &chains[ch_idx], out);
        }
        if (ch_idx%100 == 0) fprintf(stderr, "%d read chains\n", ch_idx);
    }
    ch_normalize(out, out);
}

// takes many chains of 1's and 0's (not frequencies), a position and a property to check the conditioned frequencies given that property in that position
// the property to check for conditioning must be entered as the index in the PROPERTIES list
void correlation_for_position_and_prop(Chain* chains, int n_chains, int position_index, int prop_index, Chain* out) {
    memset(out, 0, sizeof(Chain));
    for (int ch_idx = 0; ch_idx < n_chains; ch_idx++) {
        if (chains[ch_idx].aas[position_index].properties[prop_index] >= 1.-EPSILON) {
            chain_direct_sum(out, &chains[ch_idx], out);
        }
        if (ch_idx%100 == 0) fprintf(stderr, "%d read chains\n", ch_idx);
    }
    ch_normalize(out, out);
}

// applies correlations for aas and props every position and every fixed element
// **out_aacid must be CHAINLEN x N_AACIDS long
// **out_props must be CHAINLEN x N_PROPERTIES long
void first_order_correlations(Chain* input_chains, int n_chains, Chain out_aacid[][N_AACIDS], Chain out_prop[][N_PROPERTIES]) {
    for (int ch_pos = 0; ch_pos < CHAINLEN; ch_pos++) {
        for (int i = 0; i < N_AACIDS; i++) {
            correlation_for_position_and_aacid(input_chains, n_chains, ch_pos, i, &out_aacid[ch_pos][i]);
        }
        for (int i = 0; i < N_PROPERTIES; i++) {
            correlation_for_position_and_prop(input_chains, n_chains, ch_pos, i, &out_prop[ch_pos][i]);
        }
    }
}

void print_chain_matrix_to_file(Chain chs_aa[][N_AACIDS], Chain chs_pr[][N_PROPERTIES], char* filename) {
    // get FILE pointer
    FILE *f = get_file(filename, "w");

    // print header
    fprintf(f, "idxX\tidxY\t");
    for (int i = 0; i < N_AACIDS; i++) {
        fprintf(f, "%c\t", AMINOACIDS[i]);
    }
    fprintf(f, "-(PROP)\tHYDROPHOBIC\tAROMATIC\tALIPHATIC\tPOLAR\tSMALL\tMINUSCULE\tCHARGEDPLUS\tCHARGEDMINUS\n");

    for (int fdim = 0; fdim < CHAINLEN; fdim++) {
        // cycle through aminoacid-conditioned chains
        for (int sdim = 0; sdim < N_AACIDS; sdim++) {
            for (int i = 0; i < CHAINLEN; i++) {
                fprintf(f,"%d\t%d\t", fdim+1, sdim+1);
                for (int j = 0; j < N_AACIDS; j++) {
                    fprintf(f, "%.10lf\t", chs_aa[fdim][sdim].aas[i].elements[j]);
                }
                for (int j = 0; j < N_PROPERTIES; j++) {
                    fprintf(f, "%.10lf\t", chs_aa[fdim][sdim].aas[i].properties[j]);
                }
                fprintf(f, "\n");
            }
        }

        // cycle through property-conditioned chains
        for (int sdim = 0; sdim < N_PROPERTIES; sdim++) {
            for (int i = 0; i < CHAINLEN; i++) {
                fprintf(f,"%d\t%d\t", fdim+1, sdim+1);
                for (int j = 0; j < N_AACIDS; j++) {
                    fprintf(f, "%.10lf\t", chs_pr[fdim][sdim].aas[i].elements[j]);
                }
                for (int j = 0; j < N_PROPERTIES; j++) {
                    fprintf(f, "%.10lf\t", chs_pr[fdim][sdim].aas[i].properties[j]);
                }
                fprintf(f, "\n");
            }
        }
    }
    fclose(f);
}