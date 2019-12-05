#include "FrobeniusForm.h"

FrobeniusForm::FrobeniusForm(const Matrix& exMatrix) : fSize(exMatrix.size()),
                                                       a(exMatrix) {
  Matrix cur(fSize, true);
  for (int i = 0; i < fSize - 1; ++i) {
    m.push_back(cur);
    invertM.push_back(cur);
  }

  trace = getTrace();

  state = "InputMatrix";
}

bool FrobeniusForm::getFNormalFrom() {
  int i;
  for (i = fSize - 1; i > 0; --i) {
    if (fabs(a[i][i - 1]) > 1e-8) {

      generateMatriciesToMultiply(i - 1);
      a = invertM[i - 1] * a;
      a = a * m[i - 1];

    } else {
      return false;
    }
  }

  state = "FrobeniusForm";
  return true;
}

void FrobeniusForm::generateMatriciesToMultiply(int i) {
  for (int j = 0; j < fSize; ++j) {
    m[i][j][j] = 1;
    m[i][i][j] = -a[i + 1][j] / a[i + 1][i];

    invertM[i][j][j] = 1;
    invertM[i][i][j] = a[i + 1][j];

    if (i == j) {
      m[i][i][j] /= -a[i + 1][i];
    }
  }
}

float FrobeniusForm::getTrace() {
  float sum = 0;

  for (int i = 0; i < a.size(); ++i) {
    sum += a[i][i];
  }

  return sum;
}
