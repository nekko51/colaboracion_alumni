#include "head.h"

char *PROPERTIES[N_PROPERTIES] = {"-", HYDROPHOBIC, AROMATIC, ALIPHATIC, POLAR, SMALL, MINUSCULE, CHARGEDPLUS, CHARGEDMINUS};//Note that this "-" is necessary, or else entropies would return -infty
char AMINOACIDS[N_AACIDS + 1] = "-ACDEFGHIKLMNPQRSTVWY";
int PROPS_AA[N_AACIDS][N_PROPERTIES];
int CHAR_TO_INT_LUT[256];

int char_in_string(char X, char* str) {
    int i = 0;
    while (str[i] != '\0') {
        if (str[i] == X) return 1; 
        i++;
    }
    return 0;
}

//initializes the PROPS_AA matrix
void initialize_properties_matrix() {
    for(int i=0; i<N_AACIDS; i++) {
        //1 & 0 filling
        for(int j=0; j<N_PROPERTIES; j++) {
            PROPS_AA[i][j] = char_in_string(AMINOACIDS[i], PROPERTIES[j]);
        }
    }
}

//initializes the char_to_int lookup table
void initialize_char_to_int_LUT() {
    //we initialize all entries to -1
    for (int i=0; i<256; i++) {
        CHAR_TO_INT_LUT[i] = -1;
    }
    //valid AA's map to their indexes
    for (int i=0; i<N_AACIDS; i++) {
        CHAR_TO_INT_LUT[(unsigned char)AMINOACIDS[i]] = i;
    }
}

// scans next line in file f and outputs via out
void read_next_line(FILE *f, char* out) {
    if (fscanf(f, "%298s", out) != 1) {
        out[0] = '\0'; 
    }
}

// gets a char and outputs the index of that aa
int char_to_int(char X) {
    return CHAR_TO_INT_LUT[(unsigned char)X];
}

// gets an index and outputs the char of that aa
char int_to_char(int X) {
    if (X < 0 || X >= N_AACIDS) {
        fprintf(stderr, "Error: Index %d out of bounds in int_to_char; reutning '?' character\n", X);
        return('?');
    }
    return(AMINOACIDS[X]);
}

//sets the list of properties of an AA, GIVEN ITS INDEX
void properties_of_Aacid(int idx, Aacid *aaX) {
    // aaX->props = PROPS_AA[idx]; (shallow copy)
    for(int i=0; i<N_PROPERTIES; i++) {
        aaX->properties[i] = PROPS_AA[idx][i];
    }
}

//gets a char and outputs the corresponding vector aa
//completely deterministic (only outs 1. or 0.)
//also adds the binary flags for aa properties
void char_to_Aacid(char X, Aacid* aaX) {
    int idx = char_to_int(X);
    for (int i = 0; i < N_AACIDS; i++) {
        aaX->elements[i] = 0.;
        if (i == idx) aaX->elements[i] = 1.;
    }
    properties_of_Aacid(idx, aaX);
}

// reads next line and outputs it in chain formatting
void get_next_chain(FILE *f, Chain* out) {
    char line[CHAINLEN + 1]; read_next_line(f, line);
    for (int i = 0; i < CHAINLEN; i++) {
        char_to_Aacid(line[i], &out->aas[i]);
    }
}

// gives list of chains in a file, output must be at least n_lines + starting_idx long
// the first element modified is the starting_idx one
// if first file to append, starting_idx should be 0
void append_file_to_chain_vector(char* filename, int n_lines, Chain *output, int starting_idx) {
    FILE *f = get_file(filename, "r");
    for (int i = starting_idx; i < n_lines + starting_idx; i++) {
        get_next_chain(f, &output[i]);
    }
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
