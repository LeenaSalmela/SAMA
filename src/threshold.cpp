#include "threshold.h"

#include <iostream>
#include <fstream>

Thresholds::Thresholds()
{
  max_a = 500;
  th.resize(max_a+1);
  for(int i = 0; i <= max_a; i++) {
    th[i] = i+10000;
  }
}

int Thresholds::get(int a)
{
  if (a > max_a) {
    return a+10000;
  }

  return th[a];
}

void Thresholds::readThresholds(std::string filename)
{
  std::ifstream file(filename);
  std::string cell;
  
  while(!file.eof())
  {
    if (std::getline(file, cell, ' '))
    {	
      int a = std::stoi(cell);
      std::getline(file, cell, ' ');
      int t = std::stoi(cell);
      std::getline(file, cell, '\n');
      while (a > max_a) {
	th.resize(max_a*2+1);
	for(int i = max_a+1; i <= max_a*2; i++) {
	  th[i] = i+1;
	}
	max_a = 2*max_a;
      }
      th[a] = t;
    }
  }
  file.close();
}

