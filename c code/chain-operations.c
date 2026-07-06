#include "head.h"

#ifndef EPSILON 
#define EPSILON 1e-10
#endif

// sums frequencies of two aa
Aacid aacid_direct_sum(Aacid a, Aacid b) {
    Aacid out;
    for (int i = 0; i < N_AACIDS; i++) {
        out.elmts[i] = a.elmts[i] + b.elmts[i];
    }
    for (int i = 0; i < N_PROPERTIES; i++) {
        out.props[i] = a.props[i] + b.props[i];
    }
    return out;
}

// scales frequencies of an aa
Aacid aa_scale_only_aacids(Aacid aa, double scalar) {
    Aacid out;
    for (int i = 0; i < N_AACIDS; i++) {
        out.elmts[i] = aa.elmts[i] * scalar;
    }
    return out;
}

Aacid aa_scale_only_props(Aacid aa, double scalar) {
    Aacid out;
    for (int i = 0; i < N_PROPERTIES; i++) {
        out.props[i] = aa.props[i] * scalar;
    }
    return out;
}

Aacid aa_normalize_props(Aacid aa) {
    double sum = 0.;
    for (int i = 0; i < N_PROPERTIES; i++) {
        sum += aa.props[i];
    }
    if (sum > EPSILON)  return aa_scale_only_props(aa, 1/sum);
    else                return aa_scale_only_props(aa, 0.);
}

Aacid aa_normalize_aacids(Aacid aa) {
    double sum = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sum += aa.elmts[i];
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

// scales frequencies of a chain
Chain ch_normalize(Chain c) {
    Chain out;
    for (int i = 0; i < CHAINLEN; i++) {
        out.aas[i] = aa_normalize_aacids(c.aas[i]);
        out.aas[i] = aa_normalize_props(c.aas[i]);
    }
    return out;
}

// returns the schrödinger chain of a file, needs file and number of lines
Chain file_megaAacids(char *filename, int n_lines) {
    FILE *f = fopen(filename, "r");
    if (f == NULL) {printf("\nCould not open %s\n", filename);}
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
        sume += -aa.elmts[i] * log(aa.elmts[i]);
    }
    for (int i = 0; i < N_PROPERTIES; i++) {
        sump += -aa.props[i] * log(aa.props[i]);
    }
    return (Vec2){ .x = sume, .y = sump };
}

Vec2 aa_linear_entropy(Aacid aa) {
    double sume = 0., sump = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sume += aa.elmts[i] * (1 - aa.elmts[i]);
    }
    for (int i = 0; i < N_PROPERTIES; i++) {
        sump += aa.props[i] * (1 - aa.props[i]);
    }
    return (Vec2){ .x = sume, .y = sump };
}

Vec2 aa_renyi_entropy(Aacid aa, double q) {
    double sume = 0., sump = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sume += pow(aa.elmts[i],q);
    }
    for (int i = 0; i < N_PROPERTIES; i++) {
        sump += pow(aa.props[i],q);
    }
    return (Vec2){ .x = log(sume)/(1-q), .y = log(sump)/(1-q) };
}

Vec2 aa_tsallis_entropy(Aacid aa, double q) {
    double sume = 0., sump = 0.;
    for (int i = 0; i < N_AACIDS; i++) {
        sume += pow(aa.elmts[i],q);
    }
    for (int i = 0; i < N_AACIDS; i++) {
        sump += pow(aa.elmts[i],q);
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
            if (mega_chain.aas[aa].elmts[i] >= EPSILON) {
                plogp -= mega_chain.aas[aa].elmts[i] * log(mega_chain.aas[aa].elmts[i]);
                p1p += mega_chain.aas[aa].elmts[i] * (1 - mega_chain.aas[aa].elmts[i]);
                ppowq += pow(mega_chain.aas[aa].elmts[i], order);
            }
        }
        output[aa].saa = plogp; output[aa].laa = p1p; 
        output[aa].taa = (ppowq-1)/(1-order); output[aa].raa = log(ppowq)/(1-order);

        plogp = 0.; p1p = 0.; ppowq = 0.;
        for (int i = 0; i < N_PROPERTIES; i++) {
            if (mega_chain.aas[aa].props[i] >= EPSILON) {
                plogp -= mega_chain.aas[aa].props[i] * log(mega_chain.aas[aa].props[i]);
                p1p += mega_chain.aas[aa].props[i] * (1 - mega_chain.aas[aa].props[i]);
                ppowq += pow(mega_chain.aas[aa].props[i], order);
            }
        }
        output[aa].spp = plogp; output[aa].lpp = p1p; 
        output[aa].rpp = log(ppowq)/(1-order); output[aa].tpp = (ppowq-1)/(1-order);
    }
}

void print_chain(Chain c) {
    for (int i = 0; i < CHAINLEN; i++) {
        printf("%d:\t", i+1);
        for (int j = 0; j < N_AACIDS; j++) {
            printf("%.6g\t", c.aas[i].elmts[j]);
        }
        printf("\t");
        for(int j = 0; j < N_PROPERTIES; j++) {
            printf("%.6g\t", c.aas[i].props[j]);
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
            fprintf(f, "%.10lf\t", c.aas[i].elmts[j]);
        }
        for (int j = 0; j < N_PROPERTIES; j++) {
            fprintf(f, "%.10lf\t", c.aas[i].props[j]);
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
