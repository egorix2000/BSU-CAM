#include <iostream>
#include "Matrix.h"
#include "FrobeniusForm.h"
#include "PrintToFile.h"

const int N = 4;

int main() {
  bool isComplete;
  PrintToFile printer;

  do {
    Matrix matrix(N, false);
    FrobeniusForm frobeniusForm(matrix);
    printer(frobeniusForm);

    isComplete = frobeniusForm.getFNormalFrom();
    printer(frobeniusForm);

    std::cout << "Trace: " << frobeniusForm.trace;
  } while(!isComplete);

  return 0;
}