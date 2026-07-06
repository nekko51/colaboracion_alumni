#include "head.h"

char *PROPERTIES[N_PROPERTIES] = {"-", HYDROPHOBIC, AROMATIC, ALIPHATIC, POLAR, SMALL, MINUSCULE, CHARGEDPLUS, CHARGEDMINUS};//Note that this "-" is necessary, or else entropies would return -infty
char AMINOACIDS[N_AACIDS + 1] = "-ACDEFGHIKLMNPQRSTVWY";

// const int PROPS_AA[N_AACIDS][N_PROPERTIES] = {
//     {1, 0, 0, 0, 0, 0, 0, 0, 0},//DASH
//     {0, 1, 0, 0, 0, 1, 1, 0, 0},//A
//     {0, 1, 0, 0, 0, 1, 0, 0, 0},//C
//     {0, 0, 0, 0, 1, 1, 0, 0, 1},//D
//     {0, 0, 0, 0, 1, 0, 0, 0, 1},//E
//     {0, 1, 1, 0, 0, 0, 0, 0, 0},//F
//     {0, 1, 0, 0, 0, 1, 1, 0, 0},//G
//     {0, 1, 1, 0, 1, 0, 0, 1, 0},//H
//     {0, 1, 0, 1, 0, 0, 0, 0, 0},//I
//     {0, 1, 0, 0, 1, 0, 0, 1, 0},//K
//     {0, 1, 0, 1, 0, 0, 0, 0, 0},//L
//     {0, 1, 0, 0, 0, 0, 0, 0, 0},//M
//     {0, 0, 0, 0, 1, 1, 0, 0, 0},//N
//     {0, 0, 0, 0, 0, 1, 0, 0, 0},//P
//     {0, 0, 0, 0, 1, 0, 0, 0, 0},//Q
//     {0, 0, 0, 0, 1, 0, 0, 1, 0},//R
//     {0, 0, 0, 0, 1, 1, 1, 0, 0},//S
//     {0, 1, 0, 0, 1, 1, 0, 0, 0},//T
//     {0, 1, 0, 1, 0, 1, 0, 0, 0},//V
//     {0, 1, 1, 0, 1, 0, 0, 0, 0},//W
//     {0, 1, 1, 0, 1, 0, 0, 0, 0},//Y
//     {0, 0, 0, 0, 0, 0, 0, 0, 0},//UNKNOWN
// };

//initializes the PROPS_AA matrix
void initialize_props_matrix() {
    //malloc
    PROPS_AA = malloc(N_AACIDS * sizeof(int*));
    if(PROPS_AA == NULL) {
        printf("Fatal error, ran out of memory when creating PROPS_AA matrix");
        return;
    }
    
    for(int i=0; i<N_AACIDS; i++) {
        PROPS_AA[i] = malloc(N_PROPERTIES * sizeof(int));
        //1 & 0 filling
        for(int j=0; j<N_PROPERTIES; j++) {
            PROPS_AA[i][j] = char_in_string(AMINOACIDS[i], PROPERTIES[j]);
        }
    }
}

void free_props_matrix() {
    for(int i=0; i<N_AACIDS; i++) {
        free(PROPS_AA[i]);
    }
    free(PROPS_AA);
}

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

//sets the list of properties of an AA, GIVEN ITS INDEX
void properties_of_Aacid(int idx, Aacid *aaX) {
    // aaX->props = PROPS_AA[idx]; (shallow copy)
    for(int i=0; i<N_PROPERTIES; i++) {
        aaX->props[i] = PROPS_AA[idx][i];
    }
}

//gets a char and outputs the corresponding vector aa
//completely deterministic (only outs 1. or 0.)
//also adds the binary flags for aa properties
Aacid char_to_Aacid(char X) {
    Aacid aaX; int idx = char_to_int(X);
    for (int i = 0; i < N_AACIDS; i++) {
        aaX.elmts[i] = 0.;
        if (i == idx) aaX.elmts[i] = 1.;
    }
    properties_of_Aacid(idx, &aaX);

    return aaX;
}

// reads next line and outputs it in chain formatting
Chain get_next_chain(FILE *f) {
    char line[298]; read_next_line(f, line);
    Chain out;
    for (int i = 0; i < CHAINLEN; i++) {
        out.aas[i] = char_to_Aacid(line[i]);
    }
    return out;
}

// gives list of chains in a file, output must be at least n_lines + starting_idx long
// the first element modified is the starting_idx one
// if first file to append, starting_idx should be 0
void append_file_to_chain_vector(char* filename, int n_lines, Chain *output, int starting_idx) {
    FILE *f = get_file(filename, "r");
    for (int i = starting_idx; i < n_lines + starting_idx; i++) {
        output[i] = get_next_chain(f);
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
