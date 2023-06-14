//
// Created by ntdat on 6/13/2023.
//

#ifndef LANCUOIPLSSS_MINUTIAE_H
#define LANCUOIPLSSS_MINUTIAE_H
struct Minutiae {
    int x;
    int y;
    double dir;
    int type
};
struct Point{
    // x and y is coordinates, I is gray-level
    // type: 0 is not thing special, 1 is termination, 2 is bifurcation
    int i;
    int j;
    int I;
    int type;
    double dir;
};
#endif //LANCUOIPLSSS_MINUTIAE_H
