#include <iostream>
#include <math.h>

#include "LinearEquationsSystem.h"

using namespace std;

void task1() {
    LinearEquationsSystem system(12, 5, 1);
    system.fillRandomly();

    vector<vector<float>> matrix = system.A;

    system.printMatrix(4);
    cout << endl;
    cout << endl;

    system.solveWithTridiagonalMatrixAlgorithm();
    system.printYApprox();
    cout << endl;
    system.printYReal();
    cout << endl;
    cout << "Error = " << system.calculateError() << endl;


}


int main() {
    cout << "Task 1" << endl << endl;
    task1();
    return 0;
}