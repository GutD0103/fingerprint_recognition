//
// Created by ntdat on 6/13/2023.
//

#ifndef LANCUOIPLSSS_FUNCTION_H
#define LANCUOIPLSSS_FUNCTION_H
#include <math.h>
#include "Minutiae.h"
#include <stdio.h>
#include <stdlib.h>

double TangentDir(int** I, int it, int jt, int anpha);
struct Point localMax(int **I, int it, int jt, double phit, int sigma);
struct Point ridgeNearest(int **I, int _is, int _js, int sigma);
struct Point RidgeFollowing(int **I, int **T, int ic, int jc, int sigma, double muy);
int StopCriteria(int **T, int **I, int In, int inValue, int jnValue, double phin);
void minutiae(int **I);
#endif //LANCUOIPLSSS_FUNCTION_H
