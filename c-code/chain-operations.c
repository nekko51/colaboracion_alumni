#include "head.h"

// sums frequencies of two aa
Aacid aacid_direct_sum(Aacid a, Aacid b) {
    Aacid out;
    for (int i = 0; i < N_AACIDS; i++) {
        out.elements[i] = a.elements[i] + b.elements[i];
    }
    for (int i = 0; i < N_PROPERTIES; i++) {
        out.properties[i] = a.properties[i] + b.properties[i];
    }
    return out;
}

// scales frequencies of an aa
Aacid aa_scale_only_aacids(Aacid aa, double scalar) {
    Aacid out;
    for (int i = 0; i < N_AACIDS; i++) {
        out.elements[i] = aa.elements[i] * scalar;
    }
    return out;
}

Aacid aa_scale_only_properties(Aacid aa, double scalar) {
    Aacid out;
    for (int i = 0; i < N_PROPERTIES; i++) {
        out.properties[i] = aa.properties[i] * scalar;
    }
    return out;
}

Aacid aa_normalize_properties(Aacid aa) {
    double sum = 0.;
    for (int i = 0; i < N_PROPERTIES; i++) {
        sum += aa.properties[i];
    }
    if (sum > EPSILON)  return aa_scale_only_properties(aa, 1/sum);
    else                return aa_scale_only_properties(aa, 0.);
}

Aacid aa_normalize_aacids(Aacid aa) {
    double sum = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sum += aa.elements[i];
    }
    if (sum > EPSILON)  return aa_scale_only_aacids(aa, 1/sum);
    else                return aa_scale_only_aacids(aa, 0.);
}

// sums frequencies of two chains
Chain chain_direct_sum(Chain a, Chain b) {
    Chain out;
    for (int i = 0; i < CHAINLEN; i++) {
        out.aas[i] = aacid_direct_sum(a.aas[i], b.aas[i]);
    }
    return out;
}

//deep copies chain a data to chain b (this only needs to be used when both are declared using malloc? can't find good info on this online)
void chain_deep_copy(const Chain a, Chain* b) {
    for(int i=0; i<CHAINLEN; i++) {
        for(int j=0; j<N_AACIDS; j++) b->aas[i].elements[j] = a.aas[i].elements[j];
        for(int j=0; j<N_PROPERTIES; j++) b->aas[i].properties[j] = a.aas[i].properties[j];
    }
}

// scales frequencies of a chain
Chain ch_normalize(Chain c) {
    Chain out;
    for (int i = 0; i < CHAINLEN; i++) {
        out.aas[i] = aa_normalize_aacids(c.aas[i]);
        out.aas[i] = aa_normalize_properties(c.aas[i]);
    }
    return out;
}

// returns the schrödinger chain of a file, needs file and number of lines
Chain file_megaAacids(char *filename, int n_lines) {
    FILE *f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "Error: Could not open %s; returning -1.0 chain...\n", filename);
        Chain out;
        for(int i=0; i<CHAINLEN; i++) {
            for(int j=0; j<N_AACIDS; j++) out.aas[i].elements[j] = -1.;
            for(int j=0; j<N_PROPERTIES; j++) out.aas[i].properties[j] = -1.;
        }
        return(out);
    }
    Chain out = get_next_chain(f);
    for (int i = 1; i < n_lines; i++) {
        out = chain_direct_sum(out, get_next_chain(f));
        if (i%100 == 0) printf("%d read lines\n", i);
    }
    fclose(f);
    return ch_normalize(out);
}


/*entropies*/

Vec2 aa_shannon_entropy(Aacid aa) {
    double sume = 0., sump = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sume += -aa.elements[i] * log2(aa.elements[i]);
    }
    for (int i = 0; i < N_PROPERTIES; i++) {
        sump += -aa.properties[i] * log2(aa.properties[i]);
    }
    return (Vec2){ .x = sume, .y = sump };
}

Vec2 aa_linear_entropy(Aacid aa) {
    double sume = 0., sump = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sume += aa.elements[i] * (1 - aa.elements[i]);
    }
    for (int i = 0; i < N_PROPERTIES; i++) {
        sump += aa.properties[i] * (1 - aa.properties[i]);
    }
    return (Vec2){ .x = sume, .y = sump };
}

Vec2 aa_renyi_entropy(Aacid aa, double q) {
    double sume = 0., sump = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sume += pow(aa.elements[i],q);
    }
    for (int i = 0; i < N_PROPERTIES; i++) {
        sump += pow(aa.properties[i],q);
    }
    return (Vec2){ .x = log(sume)/(1-q), .y = log(sump)/(1-q) };
}

Vec2 aa_tsallis_entropy(Aacid aa, double q) {
    double sume = 0., sump = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sume += pow(aa.elements[i],q);
    }
    for (int i = 0; i < N_AACIDS; i++) {
        sump += pow(aa.elements[i],q);
    }
    return (Vec2){ .x = (1-sume)/(q-1), .y = (1-sump)/(q-1) };
}

/* calcs entropy for each entry in a chain, output must be at least CHAINLEN long
for argument ''type'':
- 'l':    linear entropy
- 'r':    Renyi entropy*
- 't':    Tsallis entropy*
- none of the above:  Shannon entropy (default)

[*] order \in [0,1), only used for these [*], otherwise ignored

output's elements are 2-vecs: (aa entropy, property entropy)
*/
void entropy_vector(Chain mega_chain, Vec2 *output, char type, double order) {
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

/*calcs all kinds of entropies for each entry in a chain,
output must be at least CHAINLEN long
order \in [0,1)
*/
void all_entropies(Chain mega_chain, Entropies *output, double order) {
    double plogp, p1p, ppowq;
    for (int aa = 0; aa < CHAINLEN; aa++) {
        plogp = 0.; p1p = 0.; ppowq = 0.;

        for (int i = 0; i < N_AACIDS; i++) {
            if (mega_chain.aas[aa].elements[i] >= EPSILON) {
                plogp -= mega_chain.aas[aa].elements[i] * log2(mega_chain.aas[aa].elements[i]);
                p1p += mega_chain.aas[aa].elements[i] * (1 - mega_chain.aas[aa].elements[i]);
                ppowq += pow(mega_chain.aas[aa].elements[i], order);
            }
        }
        output[aa].saa = plogp; output[aa].laa = p1p; 
        output[aa].taa = (ppowq-1)/(1-order); output[aa].raa = log(ppowq)/(1-order);

        plogp = 0.; p1p = 0.; ppowq = 0.;
        for (int i = 0; i < N_PROPERTIES; i++) {
            if (mega_chain.aas[aa].properties[i] >= EPSILON) {
                plogp -= mega_chain.aas[aa].properties[i] * log2(mega_chain.aas[aa].properties[i]);
                p1p += mega_chain.aas[aa].properties[i] * (1 - mega_chain.aas[aa].properties[i]);
                ppowq += pow(mega_chain.aas[aa].properties[i], order);
            }
        }
        output[aa].spp = plogp; output[aa].lpp = p1p; 
        output[aa].rpp = log(ppowq)/(1-order); output[aa].tpp = (ppowq-1)/(1-order);
    }
}

/*
should weigh every entropy type to produce a single double, weighs:
saa, laa, raa, taa
spp, lpp, rpp, tpp*/
void weigh_entropies(Entropies* input, double* output, double weighs[8]) {
    for(int i=0; i<CHAINLEN; i++) {
        // for(int j=0; j<8; j++) {
        //     output[i] = input[j].saa;
        // }
        output[i] = 0.5*input[i].saa + 0.5*input[i].spp;
    }
    return;
}

void print_chain(Chain c) {
    for (int i = 0; i < CHAINLEN; i++) {
        printf("%d:\t", i+1);
        for (int j = 0; j < N_AACIDS; j++) {
            printf("%.6g\t", c.aas[i].elements[j]);
        }
        printf("\t");
        for(int j = 0; j < N_PROPERTIES; j++) {
            printf("%.6g\t", c.aas[i].properties[j]);
        }
        printf("\n");
    }
}

void print_chain_to_file(Chain c, char* filename) {
    FILE *f = get_file(filename, "w");

    fprintf(f, "idx\t");
    for (int i = 0; i < N_AACIDS; i++) {
        fprintf(f, "%c\t", AMINOACIDS[i]);
    }
    fprintf(f, "-(PROP)\tHYDROPHOBIC\tAROMATIC\tALIPHATIC\tPOLAR\tSMALL\tMINUSCULE\tCHARGEDPLUS\tCHARGEDMINUS\n");
    for (int i = 0; i < CHAINLEN; i++) {
        fprintf(f, "%d\t", i+1);
        for (int j = 0; j < N_AACIDS; j++) {
            fprintf(f, "%.10lf\t", c.aas[i].elements[j]);
        }
        for (int j = 0; j < N_PROPERTIES; j++) {
            fprintf(f, "%.10lf\t", c.aas[i].properties[j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void print_entropies(Entropies *S) {
    for (int i = 0; i < CHAINLEN; i++) {
        printf("%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
            i+1, S[i].saa, S[i].spp, S[i].laa, S[i].lpp, S[i].raa, S[i].rpp, S[i].taa, S[i].tpp
        );
    }
}

    // Shannon AA, linear AA, Renyi AA, Tsallis AA, 
    // Shannon Props, linear Props, Renyi Props, Tsallis Props

void print_entropies_to_file(Entropies *S, char* filename) {
    FILE *f = get_file(filename, "w");
    fprintf(f, "idx\tshannonAA\tshannonPROP\tlinearAA\tlinearPROP\trenyiAA\trenyiPROP\ttsallisAA\ttsallisPROP\n");
    for (int i = 0; i < CHAINLEN; i++) {
        fprintf(f, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
            i+1, S[i].saa, S[i].spp, S[i].laa, S[i].lpp, S[i].raa, S[i].rpp, S[i].taa, S[i].tpp
        );
    }
    fclose(f);
}

/**
    @brief outputs the list of positions that are non-CDRs
    @param S the list of entropies for the chain
    @param threshold limit to which CDRs are discriminated
    @param entropy_type chosen entropy to use for the discrimination
        - 'l':    linear entropy
        - 'r':    Renyi entropy
        - 't':    Tsallis entropy
        - none of the above:  Shannon entropy (default)
    @return out_indexes positions of the non-CDRs (inputted array must be at least CHAINLEN long)
    @return n_indexes number of positions outputted
**/
void noncdr_candidates(Entropies *S, double threshold, char entropy_type, int *out_indexes, int *n_indexes) {
    double entropy; *n_indexes = 0;
    for (int pos = 0; pos < CHAINLEN; pos++) {
        if (entropy_type == 'l') {
            entropy = S[pos].laa;
        } else if (entropy_type == 'r') {
            entropy = S[pos].raa;
        } else if (entropy_type == 't') {
            entropy = S[pos].taa;
        } else {
            entropy  = S[pos].saa;
        }

        if (entropy < threshold) {
            out_indexes[*n_indexes] = pos;
            (*n_indexes)++;
        }
    }
}
