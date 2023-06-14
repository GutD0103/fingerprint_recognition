//
// Created by ntdat on 6/14/2023.
//

#ifndef LANCUOIPLSSS_ENHANCE_H
#define LANCUOIPLSSS_ENHANCE_H

#include <stdio.h>
#include <math.h>
#include "global.h"
double calculateMean(int** arr, int rows, int cols);
double calculateVariance(int** arr, int rows, int cols);
int** Normalization(int **I, int rows, int cols, int M0, int V0);

struct MaskGabor{
    int ** mask;
    int size;
};
#endif //LANCUOIPLSSS_ENHANCE_H
