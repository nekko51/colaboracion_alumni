#include "head.h"

double scalar_product(double *a, double *b, int dims) {
    double sum = 0.;
    for (int i = 0; i < dims; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

// takes many chains of 1's and 0's (not frequencies), a position and an aa to check the conditioned frequencies given that aa in that position
// the aminoacid to check for conditioning must be entered as the index in the AMINOACIDS list
Chain correlation_for_position_and_aacid(Chain* chains, int n_chains, int position_index, int aa_index) {
    Chain out = chains[0];
    for (int ch_idx = 1; ch_idx < n_chains; ch_idx++) {
        if (chains[ch_idx].aas[position_index].elements[aa_index] >= 1.-EPSILON) {
            out = chain_direct_sum(out, chains[ch_idx]);
        }
        if (ch_idx%100 == 0) printf("%d read chains\n", ch_idx);
    }
    return ch_normalize(out);
}

// takes many chains of 1's and 0's (not frequencies), a position and a property to check the conditioned frequencies given that property in that position
// the property to check for conditioning must be entered as the index in the PROPERTIES list
Chain correlation_for_position_and_prop(Chain* chains, int n_chains, int position_index, int prop_index) {
    Chain out = chains[0];
    for (int ch_idx = 1; ch_idx < n_chains; ch_idx++) {
        if (chains[ch_idx].aas[position_index].properties[prop_index] >= 1.-EPSILON) {
            out = chain_direct_sum(out, chains[ch_idx]);
        }
        if (ch_idx%100 == 0) printf("%d read chains\n", ch_idx);
    }
    return ch_normalize(out);
}

// applies correlations for aas and props every position and every fixed element
// **out_aacid must be CHAINLEN x N_AACIDS long
// **out_props must be CHAINLEN x N_PROPERTIES long
void first_order_correlations(Chain* input_chains, int n_chains, Chain out_aacid[][N_AACIDS], Chain out_prop[][N_PROPERTIES]) {
    for (int ch_pos = 0; ch_pos < CHAINLEN; ch_pos++) {
        for (int i = 0; i < N_AACIDS; i++) {
            out_aacid[ch_pos][i] = correlation_for_position_and_aacid(input_chains, n_chains, ch_pos, i);
        }
        for (int i = 0; i < N_PROPERTIES; i++) {
            out_prop[ch_pos][i] = correlation_for_position_and_prop(input_chains, n_chains, ch_pos, i);
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