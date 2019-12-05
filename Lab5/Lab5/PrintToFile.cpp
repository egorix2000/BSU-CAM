#include "FrobeniusForm.h"
#include "PrintToFile.h"

void PrintToFile::operator()(FrobeniusForm& form, int iteration) {
  int pr, w;
  std::string nameOfFile;
  nameOfFile = form.state + ".txt";
  std::ofstream out(nameOfFile);

  if (form.state == "InputMatrix") {
    pr = 0;
    w = 4;

    for (int i = 0; i < form.fSize; ++i) {
      for (int j = 0; j < form.fSize - 1; ++j) {
        out << std::setw(w) << form.a[i][j] << " &";
      }

      out << std::setw(w) << form.a[i][form.fSize - 1] << " \\\\\n";
    }
  } else if (form.state == "FrobeniusForm") {
    pr = 2;
    w = 12;

    for (int i = 0; i < form.fSize; ++i) {
      for (int j = 0; j < form.fSize - 1; ++j) {
        out << std::fixed << std::setw(w) << std::setprecision(pr)
            << form.a[i][j] << " &";
      }

      out << std::fixed << std::setw(w) << std::setprecision(pr)
          << form.a[i][form.fSize - 1]
          << " \\\\\n";
    }

    printAdditionalMatricies(form);
  }

  out.flush();
}

void PrintToFile::printAdditionalMatricies(FrobeniusForm& form) {
  for (int k = 0; k < form.fSize - 1; ++k) {
    int pr = 6, w = 11;
    std::string nameOfFile;
    nameOfFile = "M" + std::to_string(k) + ".txt";
    std::ofstream out(nameOfFile);

    for (int i = 0; i < form.fSize; ++i) {
      for (int j = 0; j < form.fSize - 1; ++j) {
        out << std::fixed << std::setw(w)
            << std::setprecision(pr) << form.m[k][i][j] << " &";
      }

      out << std::fixed << std::setw(w) << std::setprecision(pr)
          << form.m[k][i][form.fSize - 1] << " \\\\\n";
    }

    out.flush();
  }
}