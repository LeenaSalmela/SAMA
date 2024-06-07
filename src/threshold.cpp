#include "threshold.h"

#include <iostream>
#include <fstream>

Thresholds::Thresholds()
{
  max_a = 0;
}

int Thresholds::get(int a)
{
  if (a > max_a) {
    return a+1;
  }

  return th[a];
}

void Thresholds::readThresholds(std::string filename)
{
  std::ifstream file(filename);
  std::string cell;
  
  while(!file.eof())
  {
    std::getline(file, cell, ' ');
    int a = std::stoi(cell);
    std::getline(file, cell, ' ');
    int t = std::stoi(cell);
    std::getline(file, cell, '\n');
    th[a] = t;
    if (a > max_a) {
      for(int i = max_a+1; i < a; i++) {
	th[a] = a+1;
      }
    }
  }
  file.close();
}

