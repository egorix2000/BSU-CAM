#include <iostream>
#include <math.h>

#include "LinearEquationsSystem.h"

using namespace std;

void task1() {
    cout << "For k = 0" << endl;
    cout << "-------" << endl;
    LinearEquationsSystem system(11, 4, 0);
    system.fillDiagonallyRandomly();

    vector<vector<float>> matrix = system.A;

    system.printMatrix(4);
    cout << endl;
    system.printB();
    cout << endl;
    cout << endl;

    system.solveWithLdLt();
    system.printXApprox();
    cout << endl;
    system.printXReal();
    cout << endl;
    cout << "Error = " << system.calculateError() << endl;

    cout << endl << "For k = 1" << endl;
    cout << "-------" << endl;
    system.resetMatrix(matrix);
    system.changeK(1);
    system.printMatrix(4);
    cout << endl;

    system.solveWithLdLt();
    system.printXApprox();
    cout << endl;
    system.printXReal();
    cout << endl;
    cout << "Error = " << system.calculateError() << endl;

}


int main() {
    cout << "Task 1" << endl << endl;
    task1();
    return 0;
}