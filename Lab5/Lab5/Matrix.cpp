//
// Created by acer on 18.10.2019.
//

#include "Matrix.h"

Matrix::Matrix(const Matrix& matrix) : mSize(matrix.mSize) {
  a = new float* [mSize];

  for (int i = 0; i < mSize; ++i) {
    a[i] = new float[mSize];

    for (int j = 0; j < mSize; ++j) {
      a[i][j] = matrix[i][j];
    }
  }
}

Matrix::Matrix(int size, bool isNull) : mSize(size) {
  a = new float* [mSize];

  for (int i = 0; i < mSize; ++i) {
    a[i] = new float[mSize];
  }

  generateCoefs(isNull);
}

Matrix::~Matrix() {
  for (int i = 0; i < mSize; ++i) {
    delete[] a[i];
  }
}

float* Matrix::operator[](int index) const {
  return a[index];
}

int Matrix::size() const {
  return mSize;
}

Matrix& Matrix::operator* (const Matrix& other) {
  auto* result = new Matrix(mSize, true);

  for (int i = 0; i < mSize; ++i) {
    for (int j = 0; j < mSize; ++j) {
      for (int k = 0; k < mSize; ++k) {
        (*result)[i][j] += a[i][k] * other[k][j];
      }
    }
  }

  return *result;
}

Matrix& Matrix::operator= (const Matrix& other) {
  if (this == &other) return *this;

  for (int i = 0; i < mSize; ++i) {
    for (int j = 0; j < mSize; ++j) {
      a[i][j] = other[i][j];
    }
  }

  return *this;
}

void Matrix::generateCoefs(bool isNull) {
  if (isNull) {
    for (int i = 0; i < mSize; ++i) {
      for (int j = 0; j < mSize; ++j) {
        a[i][j] = 0;
      }
    }
  } else {
    std::srand(std::time(0));

    for (int i = 0; i < mSize; ++i) {
      for (int j = 0; j < mSize; ++j) {
        a[i][j] = std::rand() % 101 - 50;
      }
    }
  }
}