#include "LinearEquationsSystem.h"

LinearEquationsSystem::LinearEquationsSystem(int n, int m, int k) {
    this->m = m;
    this->n = n;
    this->k = k;
    A.resize(n, vector<float>(n));
    b.resize(n);
    xReal.resize(n);
    xApprox.resize(n);
}

void LinearEquationsSystem::resetMatrix(vector<vector<float>> newA) {
    A = newA;

    //calculate new b
    for (int i = 0; i < n; i++) {
        b[i] = 0;
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * xReal[j];
        }
    }
}

void LinearEquationsSystem::fillDiagonallyRandomly() {
    //fill vector of solutions
    for (int i = 0; i < n; i++) {
        xReal[i] = m+i;
    }

    //fill matrix without diagonal
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            A[i][j] = (-1) * rand() % 5;
            A[j][i] = 0;
        }
    }

    //fill diagonal
    for (int i = 0; i < n; i++) {
        A[i][i] = 0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                if (i > j) {
                    A[i][i] -= A[i][j];
                } else {
                    A[i][i] -= A[j][i];
                }
            }
        }
    }
    A[0][0] += pow(10, (-k));

    //calculate b
    for (int i = 0; i < n; i++) {
        b[i] = 0;
        for (int j = 0; j < n; j++) {
            if (i > j) {
                b[i] += A[i][j] * xReal[j];
            } else {
                b[i] += A[j][i] * xReal[j];
            }
        }
    }
}

void LinearEquationsSystem::fillRandomly() {
    //fill vector of solutions
    for (int i = 0; i < n; i++) {
        xReal[i] = m+i;
    }

    //fill matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = rand() % 201 - 100;
        }
    }

    //calculate b
    for (int i = 0; i < n; i++) {
        b[i] = 0;
        for (int j = 0; j < n; j++) {
            b[i] += A[i][j] * xReal[j];
        }
    }
}

void LinearEquationsSystem::changeK(int k) {
    A[0][0] -= pow(10, ( -(this->k) ));
    this->k = k;
    A[0][0] += pow(10, (-k));

    //calculate b
    for (int i = 0; i < n; i++) {
        b[i] = 0;
        for (int j = 0; j < n; j++) {
            if (i > j) {
                b[i] += A[i][j] * xReal[j];
            } else {
                b[i] += A[j][i] * xReal[j];
            }
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

void LinearEquationsSystem::printXReal() {
    cout << "xReal = ( ";
    for (int i = 0; i < n; i++) {
        cout << xReal[i] << " ";
    }
    cout << ")";
}

void LinearEquationsSystem::printXApprox() {
    cout << "xApprox = ( ";
    for (int i = 0; i < n; i++) {
        cout << xApprox[i] << " ";
    }
    cout << ")";
}

void LinearEquationsSystem::printB() {
    cout << "b = ( ";
    for (int i = 0; i < n; i++) {
        cout << b[i] << " ";
    }
    cout << ")";
}

void LinearEquationsSystem::solveWithoutSelectingMainElement() {
    //forward elimination
    float l;
    for (int step = 0; step < n-1; step++) {
        if (step == 1) {
            this->printMatrix(12);
            cout << endl;
        }
        for (int i = step + 1; i < n; i++) {
            l = A[i][step] / A[step][step];
            b[i] -= l * b[step];
            for (int j = step; j < n; j++) {
                A[i][j] -= l * A[step][j];
            }
        }
    }

    //back substitution
    float sum;
    for (int i = n-1; i >= 0; i--) {
        sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * xApprox[j];
        }
        xApprox[i] = (b[i] - sum) / A[i][i];
    }
}

void LinearEquationsSystem::solveWithSelectingMainElement() {
    //forward elimination
    float l;
    int mainRow;
    for (int step = 0; step < n-1; step++) {
        if (step == 1) {
            cout << "Main element row = " << mainRow << endl;
            this->printMatrix(12);
            cout << endl;
        }

        //select main element
        mainRow = step;
        for (int i = step + 1; i < n; i++) {
            if (abs(A[i][step]) > abs(A[mainRow][step]))
                mainRow = i;
        }

        if (abs(A[mainRow][step]) == 0) {
            continue;
        }

        for (int i = step; i <= n; i++) {
            swap (A[mainRow][i], A[step][i]);
        }
        swap (b[mainRow], b[step]);

        for (int i = step + 1; i < n; i++) {
            l = A[i][step] / A[step][step];
            b[i] -= l * b[step];
            for (int j = step; j < n; j++) {
                A[i][j] -= l * A[step][j];
            }
        }
    }

    //back substitution
    float sum;
    for (int i = n-1; i >= 0; i--) {
        sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * xApprox[j];
        }
        xApprox[i] = (b[i] - sum) / A[i][i];
    }
}

void LinearEquationsSystem::solveWithLdLt() {
    //LdLt decomposition
    float l;
    float* t = new float[n];
    for (int step = 0; step < n-1; step++) {
        for (int i = step + 1; i < n; i++) {
            t[i] = A[i][step];
            A[i][step] /= A[step][step];

            for (int j = step + 1; j <= i; j++) {
                A[i][j] -= A[i][step] * t[j];
            }
        }
    }
    printMatrix(13);
    cout << endl;
    delete [] t;

    //solve Ly = b where y = DL'x
    float sum;
    float* yApprox = new float[n];
    for (int i = 0; i < n; i++) {
        sum = 0;
        for (int j = 0; j < i; j++) {
            sum += A[i][j] * yApprox[j];
        }
        yApprox[i] = b[i] - sum;
    }

    //solve DL'x = y
    //back substitution
    sum = 0;
    for (int i = n-1; i >= 0; i--) {
        sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][i] * A[j][i] * xApprox[j]; //D[i][i] * L'[i][j] * xApprox[j]
        }
        xApprox[i] = (yApprox[i] - sum) / A[i][i];
    }

    delete [] yApprox;
}

float LinearEquationsSystem::calculateError() {
    float error;
    float normError = 0;
    float xRealNorm = 0;

    for (int i = 0; i < n; i++) {
        if (abs(xReal[i]) > xRealNorm) {
            xRealNorm = abs(xReal[i]);
        }
    }

    for (int i = 0; i < n; i++) {
        if (abs(xReal[i] - xApprox[i]) > normError) {
            normError = abs(xReal[i] - xApprox[i]);
        }
    }
    error = normError / xRealNorm;
    return error;
}