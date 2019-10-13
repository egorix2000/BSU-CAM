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
    vector<float> b;
    vector<float> xReal;
    vector<float> xApprox;
    LinearEquationsSystem(int n, int m, int k = 0);
    void resetMatrix(vector<vector<float>> newA);
    void fillDiagonallyRandomly();
    void fillRandomly();
    void changeK(int k);
    void printMatrix(int width);
    void printXReal();
    void printXApprox();
    void printB();
    void solveWithoutSelectingMainElement();
    void solveWithSelectingMainElement();
    void solveWithLdLt();
    float calculateError();
};


#endif //CAM_LINEAREQUATIONSSYSTEM_H
