#ifndef LIB_CONV_H
#define LIB_CONV_H


#include <stdlib.h>
#include <math.h>

typedef struct {
    int num_input_bit;
    int num_output_bit;
    int num_state;
    int **next_state;
    int ***output;
} conv_coder;


// 畳み込み符号化(拘束長3)
conv_coder construct_encoder (int type) {
    static conv_coder coder;
    int i, j;

    coder.num_input_bit = 1;
    coder.num_output_bit = 2;
    coder.num_state = 4;

    coder.next_state = (int **)malloc(coder.num_state * sizeof(int *));
    coder.output = (int ***)malloc(coder.num_output_bit * sizeof(int **));

    for (i = 0; i < coder.num_output_bit; i++) {
        coder.output[i] = (int **)malloc(coder.num_state * sizeof(int *));

        for (j = 0; j < coder.num_state; j++) {
            coder.output[i][j] = (int *)malloc((int)pow(2.0, coder.num_input_bit) * sizeof(int));
            coder.next_state[j] = (int *)malloc((int)pow(2.0, coder.num_input_bit) * sizeof(int));
        }
    }

    return coder;
}


#endif
