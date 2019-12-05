//
// Created by acer on 18.10.2019.
//

#ifndef LAB4_MATRIX_H
#define LAB4_MATRIX_H
#include <ctime>
#include <cstdlib>

class Matrix {
 public:
  Matrix(const Matrix&);
  Matrix(int, bool);
  ~Matrix();

  float* operator[] (int) const;
  int size() const;

  Matrix& operator* (const Matrix&);
  Matrix& operator= (const Matrix&);

 protected:
  void generateCoefs(bool);

  int mSize;
  float** a;
};

#endif //LAB4_MATRIX_H
