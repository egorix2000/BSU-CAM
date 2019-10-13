#ifndef CAM_LINEAREQUATIONSSYSTEM_H
#define CAM_LINEAREQUATIONSSYSTEM_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>

using namespace std;

class LinearEquationsSystem {
public:
    int n;
    int m;
    int k;
    vector<vector<float>> A;
    vector<float> d;
    vector<float> yReal;
    vector<float> yApprox;
    LinearEquationsSystem(int n, int m, int k = 0);
    void fillRandomly();
    void printMatrix(int width);
    void printYReal();
    void printYApprox();
    void printD();
    void solveWithTridiagonalMatrixAlgorithm();
    float calculateError();
};


#endif //CAM_LINEAREQUATIONSSYSTEM_H
