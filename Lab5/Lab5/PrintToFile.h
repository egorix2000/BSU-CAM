#pragma once
#include <string>
#include <iomanip>
#include <fstream>

class FrobeniusForm;

struct PrintToFile {
  void operator() (FrobeniusForm&, int iteration = 0);

 private:
  void printAdditionalMatricies(FrobeniusForm&);
};