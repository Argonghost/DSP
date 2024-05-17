#include <stdio.h>
#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>

// GPS Shift Register function:  the function "shift" emulates a Linear Feedback Shift Register
int shift(int *reg, int *feedback, int feedback_len, int *output, int output_len) {
    int out = 0;
    for (int i = 0; i < output_len; i++) {
        out += reg[output[i] - 1];
    }
    out %= 2;

    int fb = 0;
    for (int i = 0; i < feedback_len; i++) {
        fb += reg[feedback[i] - 1];
    }
    fb %= 2;

    for (int i = 9; i > 0; i--) {
        reg[i] = reg[i - 1];
    }
    reg[0] = fb;

    return out;
}

// Build the CA code (PRN) for a given satellite ID
int PRN(int sv, int *ca) {
    int G1[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int G2[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    int feedback1[2] = {3, 10};
    int output1[1] = {10};

    int feedback2[6] = {2, 3, 6, 8, 9, 10};
    int output2[2];
    switch (sv) {
        case 1: output2[0] = 2; output2[1] = 6; break;
        case 2: output2[0] = 3; output2[1] = 7; break;
        case 3: output2[0] = 4; output2[1] = 8; break;
        case 4: output2[0] = 5; output2[1] = 9; break;
        case 5: output2[0] = 1; output2[1] = 9; break;
        case 6: output2[0] = 2; output2[1] = 10; break;
        case 7: output2[0] = 1; output2[1] = 8; break;
        case 8: output2[0] = 2; output2[1] = 9; break;
        case 9: output2[0] = 3; output2[1] = 10; break;
        case 10: output2[0] = 2; output2[1] = 3; break;
        case 11: output2[0] = 3; output2[1] = 4; break;
        case 12: output2[0] = 5; output2[1] = 6; break;
        case 13: output2[0] = 6; output2[1] = 7; break;
        case 14: output2[0] = 7; output2[1] = 8; break;
        case 15: output2[0] = 8; output2[1] = 9; break;
        case 16: output2[0] = 9; output2[1] = 10; break;
        case 17: output2[0] = 1; output2[1] = 4; break;
        case 18: output2[0] = 2; output2[1] = 5; break;
        case 19: output2[0] = 3; output2[1] = 6; break;
        case 20: output2[0] = 4; output2[1] = 7; break;
        case 21: output2[0] = 5; output2[1] = 8; break;
        case 22: output2[0] = 6; output2[1] = 9; break;
        case 23: output2[0] = 1; output2[1] = 3; break;
        case 24: output2[0] = 4; output2[1] = 6; break;
        case 25: output2[0] = 5; output2[1] = 7; break;
        case 26: output2[0] = 6; output2[1] = 8; break;
        case 27: output2[0] = 7; output2[1] = 9; break;
        case 28: output2[0] = 8; output2[1] = 10; break;
        case 29: output2[0] = 1; output2[1] = 6; break;
        case 30: output2[0] = 2; output2[1] = 7; break;
        case 31: output2[0] = 3; output2[1] = 8; break;
        default: output2[0] = 4; output2[1] = 9; break; // Default to satellite 32
    }

    for (int i = 0; i < 1023; i++) {
        int g1 = shift(G1, feedback1, 2, output1, 1);
        int g2 = shift(G2, feedback2, 6, output2, 2);
        ca[i] = (g1 + g2) % 2;
    }

    return 0;
}

// The below function, get_code, takes in a given prn code, chipping rate, and a given ADC sampling rate n, 
// in order to have the prn repeat itself in accordance with an ADC.

int* get_code(int* code, int n, int code_len) {
    int* gc = malloc(n * code_len * sizeof(int));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < code_len; j++) {
            gc[i*code_len + j] = code[j];
        }
    }
    return gc;
}

int main() {
    int sat_24[1023];
    PRN(1, sat_24);
    printf("PRN code for satellite :\n");
    printf("------------------------------\n");
    printf("{");
    for (int i = 0; i < 1023; i++) {
        printf("%d,", sat_24[i]);
    }

    printf("}");

    int *new = get_code(sat_24, 4, 1023);
    for(int i = 0; i < 1023*4; i++){
        printf("%d", new[i]);
    }
    free(new);
    return 0;
}