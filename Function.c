//
// Created by ntdat on 6/13/2023.
//
#include "Function.h"

// Determine the tangent of the Point (it,jt)
double TangentDir(int** I, int it, int jt, int anpha)
{
    double A = 0;
    double B = 0;
    double C = 0;
    int halfAnpha = (anpha - 1) / 2;
    for (int h = it - halfAnpha; h <= it + halfAnpha; h++)
    {
        for (int k = jt - halfAnpha; k <= jt + halfAnpha; k++)
        {
            double ahk = (-I[h - 1][k - 1] + I[h - 1][k + 1] + I[h + 1][k + 1] - I[h + 1][k - 1]) / 4;
            double bhk = (-I[h - 1][k - 1] - I[h - 1][k + 1] + I[h + 1][k + 1] + I[h + 1][k - 1]) / 4;
            A += pow(ahk, 2);
            B += pow(bhk, 2);
            C += ahk * bhk;
        }
    }

    double t1, t2;
    if (C > 0)
    {
        t1 = 1;
        t2 = (B - A) / (2 * C) - sqrt(pow((B - A) / (2 * C), 2) + 1);
    }
    else if (C < 0)
    {
        t1 = 1;
        t2 = (B - A) / (2 * C) + sqrt(pow((B - A) / (2 * C), 2) + 1);
    }
    else if (A <= B)
    {
        t1 = 1;
        t2 = 0;
    }
    else
    {
        t1 = 0;
        t2 = 1;
    }

    double phit;
    if (t1 == 0)
    {
        phit = M_PI / 2;
    }
    else
    {
        phit = atan(t2 / t1);
    }

    return phit;
}

// This function is used for find a local maximum
struct Point localMax(int **I, int it, int jt, double phit, int sigma) {
    double d[] = { 1.0 / 23, 2.0 / 23, 5.0 / 23, 7.0 / 23, 5.0 / 23, 2.0 / 23, 1.0 / 23 };
    double LM[sigma * 2 + 1];
    int io[sigma * 2 + 1];
    int jo[sigma * 2 + 1];
    struct Point result = {0,0,0,0,0};
    result.I = 0;
    result.i = -1;
    result.j = -1;

    for (int x = -sigma; x <= sigma; x++) {
        int i = (int)(it + x * sin(phit + M_PI / 2));
        int j = (int)(jt + x * cos(phit + M_PI / 2));
        int iPrev = (int)(i - sin(phit));
        int jPrev = (int)(j - cos(phit));
        int iNext = (int)(i + sin(phit));
        int jNext = (int)(j + cos(phit));

        double Io = (I[iPrev][jPrev] + I[i][j] + I[iNext][jNext]) / 3.0;
        LM[x + sigma] = Io;
        io[x + sigma] = i;
        jo[x + sigma] = j;
    }

    double LMn[sigma * 2 + 1];
    for (int x = 4; x <= 2 * sigma - 2; x++) {
        LMn[x] = LM[x - 3] * d[0] + LM[x - 2] * d[1] + LM[x - 1] * d[2] + LM[x] * d[3] +
                 LM[x + 1] * d[4] + LM[x + 2] * d[5] + LM[x + 3] * d[6];

        if (result.I < LMn[x]) {
            result.i = io[x];
            result.j = jo[x];
            result.I = (int)LMn[x];
        }
    }
    return result;
}

struct Point ridgeNearest(int **I, int _is, int _js, int sigma) {

    double phiC = TangentDir(I, _is, _js, 9);
    int i1 = round(_is - (sigma - 3) * cos(phiC + M_PI / 2));
    int j1 = round(_js - (sigma - 3) * sin(phiC + M_PI / 2));

    struct Point _p1 = localMax(I, i1, j1, phiC, sigma);

    int i2 = round(_is + (sigma - 3) * cos(phiC + M_PI / 2));
    int j2 = round(_js + (sigma - 3) * sin(phiC + M_PI / 2));

    struct Point _p2 = localMax(I, i2, j2, phiC, sigma);

    int d1 = (_p1.i - _is) * (_p1.i - _is) + (_p1.j - _js) * (_p1.j - _js);
    int d2 = (_p2.i - _is) * (_p2.i - _is) + (_p2.j - _js) * (_p2.j - _js);

    struct Point result;
    if (d1 < d2) {
        result = _p1;
    } else {
        result = _p2;
    }

    return result;
}

struct Point RidgeFollowing(int **I, int **T, int ic, int jc, int sigma, double muy) {
    int m = sizeof(I) / sizeof(I[0]);
    int n = sizeof(I[0]) / sizeof(I[0][0]);

    int io = ic;
    int jo = jc;

    double muy0 = (muy > 0) ? muy : -muy;
    double object0 = (muy > 0) ? 0 : M_PI;

    double phi00 = TangentDir(I, io, jo, 9);
    double phi0 = phi00 + object0;

    int flag_esc = 1;
    int type;

    while (flag_esc) {
        int it = (int) round(io + muy0 * sin(phi0));
        int jt = (int) round(jo + muy0 * cos(phi0));

        if (it > 26 && it < m - 25 && jt > 26 && jt < n - 25) {

            double phit = TangentDir(I, it, jt, 9) + object0;

            struct Point _localMax = localMax(I, it, jt, phit, sigma);

            if (_localMax.i == io && _localMax.j == jo) {
                flag_esc = 0;
            } else {

                double phinn = TangentDir(I, _localMax.i, _localMax.j, 9);
                double phin = phinn + object0;
                type = StopCriteria(T, I, _localMax.I, _localMax.i, _localMax.j, phin);
                io = _localMax.i;
                jo = _localMax.j;
                phi0 = phin;

                if(type != 1){
                    flag_esc = 0;
                }
            }
        } else {
            flag_esc = 0;
            type = 1;
        }
    }
    struct Point output = {0,0,0,0,0};
    output.i = io;
    output.j = jo;
    output.type = type;

    return output;
}
int StopCriteria(int **T, int **I, int In, int inValue, int jnValue, double phin) {
        double phinDegrees = phin * 180 / M_PI;
        int I_threshold = 20;
        int flag = 0; // flag for check this point is a bifurcation
        int output; // output of this function
        int a, b, c, d;
        if (((phinDegrees > -90) && (phinDegrees <= -67.5)) || ((phinDegrees > 247.5) && (phinDegrees <= 270))) {
            if (T[inValue - 1][jnValue] == 1)
                flag = 1;

            a = 6;
            b = 4;
            c = 4;
            d = 4;
        } else if ((phinDegrees > -67.5) && (phinDegrees <= -22.5)) {
            if (T[inValue - 1][jnValue + 1] == 1)
                flag = 1;

            a = 6;
            b = 4;
            c = 4;
            d = 6;
        } else if ((phinDegrees > -22.5) && (phinDegrees <= 22.5)) {
            if (T[inValue][jnValue + 1] == 1)
                flag = 1;

            a = 4;
            b = 4;
            c = 4;
            d = 6;
        } else if ((phinDegrees > 22.5) && (phinDegrees <= 67.5)) {
            if (T[inValue + 1][jnValue + 1] == 1)
                flag = 1;

            a = 4;
            b = 6;
            c = 4;
            d = 6;
        } else if ((phinDegrees > 67.5) && (phinDegrees <= 112.5)) {
            if (T[inValue + 1][jnValue] == 1)
                flag = 1;

            a = 4;
            b = 6;
            c = 4;
            d = 4;
        } else if ((phinDegrees > 112.5) && (phinDegrees <= 157.5)) {
            if (T[inValue + 1][jnValue - 1] == 1)
                flag = 1;

            a = 4;
            b = 6;
            c = 6;
            d = 4;
        } else if ((phinDegrees > 157.5) && (phinDegrees <= 202.5)) {
            if (T[inValue][jnValue - 1] == 1)
                flag = 1;

            a = 4;
            b = 4;
            c = 6;
            d = 4;
        } else {
            if (T[inValue - 1][jnValue + 1] == 1)
                flag = 1;

            a = 6;
            b = 4;
            c = 6;
            d = 4;
        }


        if (flag) {
            int label = 1;
            // If in the range around point X have 1 or more Minutiae are founded,
            // it's an error minutiae
            for (int i = inValue - a; i <= inValue + b; i++) {
                for (int j = jnValue - c; j <= jnValue + d; j++) {
                    if ((T[i][j] == 2) || (T[i][j] == 3) || (T[i][j] == -1)) {
                        T[i][j] = -1; // A error minutiae had been found
                        label = 0;
                    }
                }
            }
            // if this Minutiae is not error and these points around 20 pixel is not the ridge, so it's a bifurcation
            if (label && (T[inValue - 20][jnValue] != 1) && (T[inValue + 20][jnValue] != 1)
                && (T[inValue][jnValue - 20] != 1) && (T[inValue][jnValue + 20] != 1))
                T[inValue][jnValue] = 3; // A bifurcation

            // stop
        } else {
            for (int i = inValue - 1; i <= inValue + 1; i++) {
                for (int j = jnValue - 1; j <= jnValue + 1; j++) {
                    T[i][j] = 1;
                }
            }

            if ((In < I_threshold) && (T[inValue - 20][jnValue] != 1) && (T[inValue + 20][jnValue] != 1) &&
                (T[inValue][jnValue - 20] != 1) && (T[inValue][jnValue + 20] != 1)) {
                int label = 1;
                for (int i = inValue - a; i <= inValue + b; i++) {
                    for (int j = jnValue - c; j <= jnValue + d; j++) {
                        if ((T[i][j] == 2) || (T[i][j] == 3) || (T[i][j] == -1)) {
                            T[i][j] = -2;
                            label = 0;
                        }
                    }
                }
                if (label)
                    T[inValue][jnValue] = 2; // A termination
            }

        }
        return T[inValue][jnValue];
}
void minutiae(int **I) {

        int m = sizeof(I) / sizeof(I[0]);
        int n = sizeof(I[0]) / sizeof(I[0][0]);

        int v = 4;
        int sigma = 5;

        int **T = (int **) malloc(m * sizeof(int *));
        for (int i = 0; i < m; i++) {
            T[i] = (int *) malloc(n * sizeof(int));
        }


        for (int is = 26; is < m - 25; is += v) {
            for (int js = 26; js < n - 25; js += v) {

                    struct Point _ridgeNearest = ridgeNearest(I, is, js, sigma);

                    int label = 1;

                    for (int i = _ridgeNearest.i - 4; i <= _ridgeNearest.i + 4; i++) {
                        for (int j = _ridgeNearest.j - 4; j <= _ridgeNearest.j + 4; j++) {
                            if (T[i][j] == 1) {
                                label = 0;
                            }
                        }
                    }

                    if (label) {
                       struct Point temp1 = RidgeFollowing(I, T, _ridgeNearest.i, _ridgeNearest.j, sigma, 2.4);
                       struct Point temp2 = RidgeFollowing(I, T, _ridgeNearest.i, _ridgeNearest.j, sigma, -2.4);
                    }
            }
        }
}
