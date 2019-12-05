#pragma once
#include <string>
#include <cmath>
#include <vector>
#include "PrintToFile.h"
#include "Matrix.h"

class FrobeniusForm{
 public:

  explicit FrobeniusForm(const Matrix&);
  float trace;

  bool getFNormalFrom();

 private:
  int fSize;
  std::string state;
  Matrix a;
  std::vector<Matrix> m;
  std::vector<Matrix> invertM;

  float getTrace();
  void generateMatriciesToMultiply(int);

  friend class PrintToFile;
};