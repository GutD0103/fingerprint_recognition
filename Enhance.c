// Created by ntdat on 6/14/2023.
//


// Hàm tính trung bình của mảng 2 chiều

#include "Enhance.h"

double calculateMean(int** arr, int rows, int cols) {
    double sum = 0;
    int count = 0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            sum += arr[i][j];
            count++;
        }
    }

    return sum / count;
}

// Hàm tính phương sai của mảng 2 chiều
double calculateVariance(int** arr, int rows, int cols) {
    double mean = calculateMean(arr, rows, cols);
    double sum = 0;
    int count = 0;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double diff = arr[i][j] - mean;
            sum += diff * diff;
            count++;
        }
    }

    return sum / count;
}

int** Normalization(int **I, int rows, int cols, int M0, int V0){
    double  M = calculateMean(I,rows,cols);
    double  V = calculateVariance(I,rows,cols);
    int** G = create_2d_array(rows,cols);
    for(int i = 0; i < rows ; i++){
        for (int j = 0; j < cols; ++j) {
            if(I[i][j] > M){
                G[i][j] =  M0 + sqrt(V0 * pow((I[i][j] - M), 2) / V);
            }else{
                G[i][j] =  M0 - sqrt(V0 * pow((I[i][j] - M), 2) / V);
            }
        }
    }
    return G;
}

double GetDirectionMatrix(int widthSquare, int width, int height, int** I) {
    int i, j, x, y;
    int Ax, Ay, Axy;
    int Gxx, Gyy, Gxy;
    int Bx, By;

    double directMatrix;

            Ax = 0; Ay = 0; Axy = 0;
            Gxx = 0; Gyy = 0; Gxy = 0;
            Bx = 0; By = 0;

            for (j = y - widthSquare; j < y + widthSquare - 1; j++) {
                for (i = x - widthSquare; i < x + widthSquare - 1; i++) {
                    Bx = ((I[i + 2][j] + 2 * I[i + 2][j + 1] + I[i + 2][j + 2] - I[i][j] - 2 * I[i][j + 1] - I[i][j + 2]));
                    By = ((I[i][j + 2] + 2 * I[i + 1][j + 2] + I[i + 2][j + 2] - I[i][j] - 2 * I[i + 1][j] - I[i + 2][j]));
                    Ax += Bx * Bx;
                    Ay += By * By;
                    Axy += Bx * By;
                }
            }

            Gxx = Ax;
            Gyy = Ay;
            Gxy = Axy;

            directMatrix = M_PI / 2 - 0.5 * atan2(2 * Gxy, Gxx - Gyy);

    return directMatrix;
}
void MaskGabor(int widthSquare, double filterDirect, double f, int fi, double** mask, int* width) {
    int size = 2 * widthSquare + 1;
    int** rtMask = create_2d_array(size,size);

    for (int x = 0; x < size; x++) {
        for (int y = 0; y < size; y++) {
            double x1 = sin(filterDirect) * (x - widthSquare) + cos(filterDirect) * (y - widthSquare);
            double y1 = sin(filterDirect) * (y - widthSquare) - cos(filterDirect) * (x - widthSquare);
            rtMask[x][y] = exp(-0.5 * (pow(x1, 2) / pow(fi, 2) + pow(y1, 2) / pow(fi, 2))) * cos(2 * M_PI * f * x1);
        }
    }
    struct MaskGabor result;
    result.mask = rtMask;
    result.size = size;
}
void ToFiltring(int widthSquare, int f, int fi, int width, int height, int ** I) {
        double pointValue = 0;
        for (int x = 0; x < width - 2 * widthSquare - 1; x++) {
            for (int y = 0; y < height - 2 * widthSquare - 1; y++) {
                double** mask;
                if (FingerCompare_maskNumber > 0) {
                    mask = GetMaskFilter(direct[x][y]);
                } else {
                    mask = GetMaskFilter(direct[x][y], widthSquare, 1.0 / f, fi);
                }
                for (int i = 0; i < 2 * widthSquare + 1; i++) {
                    for (int j = 0; j < 2 * widthSquare + 1; j++) {
                        pointValue += mask[i][j] * I[i + x][j + y];
                    }
                }
                if (pointValue < 0)
                    pointValue = 0;
                if (pointValue > 255)
                    pointValue = 255;
                I[x][y] = (int)pointValue;
            }
        }
}

double** GetMaskFilter(double filterDirect) {
    MaskGabor mask;
    for (int i = 0; i < FingerCompare_maskNumber; i++) {
        if (filterDirect >= i * M_PI / FingerCompare_maskNumber && filterDirect < (i + 1) * M_PI / FingerCompare_maskNumber) {
            mask = (MaskGabor)FingerCompare_MaskGaborCollection[i];
            return mask.mask;
        }
    }
    mask = (MaskGabor)FingerCompare_MaskGaborCollection[FingerCompare_maskNumber - 1];
    return mask.mask;
}

double** GetMaskFilter(double filterDirect, int widthSquare, double f, int fi) {
    int size = 2 * widthSquare + 1;
    double** rtMask = (double**)malloc(size * sizeof(double*));
    for (int i = 0; i < size; i++) {
        rtMask[i] = (double*)malloc(size * sizeof(double));
    }

    for (int x = 0; x < size; x++) {
        for (int y = 0; y < size; y++) {
            double x1 = sin(filterDirect) * (x - widthSquare) + cos(filterDirect) * (y - widthSquare);
            double y1 = sin(filterDirect) * (y - widthSquare) - cos(filterDirect) * (x - widthSquare);
            rtMask[x][y] = exp(-0.5 * (pow(x1, 2) / pow(fi, 2) + pow(y1, 2) / pow(fi, 2))) * cos(2 * M_PI * f * x1);
        }
    }
    return rtMask;
}
