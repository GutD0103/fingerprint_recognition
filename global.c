//
// Created by ntdat on 6/14/2023.
//
#include "global.h"
int** create_2d_array(int rows, int cols) {
    int** arr = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++) {
        arr[i] = (int*)malloc(cols * sizeof(int));
    }
    return arr;
}

void destroy_2d_array(int** arr, int rows) {
    for (int i = 0; i < rows; i++) {
        free(arr[i]);
    }
    free(arr);
}