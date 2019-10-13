#include "LinearEquationsSystem.h"

LinearEquationsSystem::LinearEquationsSystem(int n, int m, int k) {
    this->m = m;
    this->n = n;
    this->k = k;
    A.resize(n, vector<float>(n));
    d.resize(n);
    yReal.resize(n);
    yApprox.resize(n);
}

void LinearEquationsSystem::fillRandomly() {
    //fill vector of solutions
    for (int i = 0; i < n; i++) {
        yReal[i] = i+1;
    }

    //fill matrix
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            A[i][i] = m;
            A[i][i+1] = m - 1;
        }
        else if (i == n-1) {
            A[i][i] = m + k + i - 1;
            A[i][i-1] = -k;
        }
        else {
            A[i][i] = m + k + i - 1;
            A[i][i-1] = -k;
            A[i][i+1] = m + i - 1;
        }
    }

    //calculate b
    for (int i = 0; i < n; i++) {
        d[i] = 0;
        for (int j = 0; j < n; j++) {
            d[i] += A[i][j] * yReal[j];
        }
    }
}

void LinearEquationsSystem::printMatrix(int width) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(width) << A[i][j];
        }
        cout << endl;
    }
}

void LinearEquationsSystem::printYReal() {
    cout << "yReal = ( ";
    for (int i = 0; i < n; i++) {
        cout << yReal[i] << " ";
    }
    cout << ")";
}

void LinearEquationsSystem::printYApprox() {
    cout << "yApprox = ( ";
    for (int i = 0; i < n; i++) {
        cout << setprecision(10) << yApprox[i] << " ";
    }
    cout << ")";
}

void LinearEquationsSystem::printD() {
    cout << "d = ( ";
    for (int i = 0; i < n; i++) {
        cout << d[i] << " ";
    }
    cout << ")";
}

void LinearEquationsSystem::solveWithTridiagonalMatrixAlgorithm() {
    //forward sweep
    float* alpha = new float[n-1];
    float* beta = new float[n];
    float denominatorTemp;

    alpha[0] = A[0][1] / A[0][0];
    beta[0] = d[0] / A[0][0];
    for (int i = 1; i < n - 1; i++) {
        denominatorTemp = A[i][i] - A[i][i-1] * alpha[i-1];
        alpha[i] = A[i][i+1] / denominatorTemp;
        beta[i] = (d[i] - A[i][i-1] * beta[i-1]) / denominatorTemp;
    }

    beta[n - 1] = (d[n-1] - A[n-1][n-2] * beta[n-2]) / (A[n-1][n-1] - A[n-1][n-2] * alpha[n-2]);

    //print
    cout << "alpha: ( ";
    for (int i = 0; i < n-1; i++) {
        cout << alpha[i] << " ";
    }
    cout << ")" << endl;

    cout << "beta: ( ";
    for (int i = 0; i < n-1; i++) {
        cout << beta[i] << " ";
    }
    cout << ")" << endl;

    //back substitution
    yApprox[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        yApprox[i] = beta[i] - alpha[i] * yApprox[i+1];
    }
}

float LinearEquationsSystem::calculateError() {
    float error;
    float normError = 0;
    float yRealNorm = 0;

    for (int i = 0; i < n; i++) {
        yRealNorm += abs(yReal[i]);
    }

    for (int i = 0; i < n; i++) {
        normError += abs(yReal[i] - yApprox[i]);
    }
    error = normError / yRealNorm;
    return error;
}