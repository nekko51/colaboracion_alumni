#include "head.h"

char *PROPERTIES[N_PROPERTIES] = {"-", HYDROPHOBIC, AROMATIC, ALIPHATIC, POLAR, SMALL, MINUSCULE, CHARGEDPLUS, CHARGEDMINUS};

// scans next line in file f and outputs via out
void read_next_line(FILE *f, char* out) {
    fscanf(f,"%298s", out);
}

// gets a char and outputs the index of that aa
int char_to_int(char X) {
    char *aa = AACIDS;
    for (int i = 0; i < N_AACIDS; i++) {
        if (X == aa[i]) return i;
    }
    return -1;
}

int char_in_string(char X, char* str) {
    int i = 0;
    while (str[i] != '\0') {
        if (str[i] == X) return 1; 
        i++;
    }
    return 0;
}

// gets a char and outputs the corresponding vector aa
//completely deterministic (only outs 1. or 0.)
//also adds the binary flags for aa properties
Aacid char_to_Aacid(char X) {
    Aacid aaX; int idx = char_to_int(X);
    for (int i = 0; i < N_AACIDS; i++) {
        aaX.elmts[i] = 0.;
        if (i == idx) aaX.elmts[i] = 1.;
    }

    for (int i = 0; i < N_PROPERTIES; i++) {
        if (char_in_string(X, PROPERTIES[i]) == 1)  aaX.props[i] = 1.;
        else                                        aaX.props[i] = 0.;
    }

    return aaX;
}

// reads next line and outputs it in chain formatting
Chain get_nex_chain(FILE *f) {
    char line[298]; read_next_line(f, line);
    Chain out;
    for (int i = 0; i < CHAINLEN; i++) {
        out.aas[i] = char_to_Aacid(line[i]);
    }
    return out;
}





// int main() {
//     FILE *f = fopen("../seqs/learn_human.txt", "r");
//     if (f == NULL) return 1;

//     char line[297];

//     for (int i = 0; i < 3; i++) {
//     read_next_line(f, line);
//     printf("%s\n", line);}

//     fclose(f);
//     return 0;
    
// }
