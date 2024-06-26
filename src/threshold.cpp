#include "threshold.h"

#include <iostream>
#include <fstream>
#include <tuple>
#include <math.h>
#include <stdlib.h>

Thresholds::Thresholds()
{
  max_a = 500;
  th.resize(max_a+1);
  for(int i = 0; i <= max_a; i++) {
    th[i] = i+10000;
  }
  prob.resize(max_a+1);
  for(int i = 0; i <= max_a; i++) {
    prob[i] = 0.0;
  }
}

int Thresholds::get(int a)
{
  if (a > max_a) {
    return a+10000;
  }

  return th[a];
}

void Thresholds::computeThresholds(std::string filename, double epsilon) {
  std::ifstream file(filename);
  std::string cell;

  std::vector<double *> repeatP;
  int max_cov = 10;

  std::cout << epsilon << std::endl;
  
  repeatP.resize(max_a+1);
  for(int a = 0; a <= max_a; a++) {
    repeatP[a] = NULL;
  }
  
  while(!file.eof())
  {
    if (std::getline(file, cell, '\t'))
    {	
      double *row = new double[5];
      int cov = std::stoi(cell);
      std::getline(file, cell, '\t');  // 0 copies
      double s = 0.0;
      for(int r = 1; r <= 5; r++) {
	std::getline(file, cell, '\t');  // r copies
	double p = std::strtod(cell.c_str(), NULL);
	row[r-1] = p;
	s+= p;
      }
      std::getline(file, cell, '\n');
      for(int r = 1; r <= 5; r++) {
	row[r-1] = row[r-1]/s;
      }
      if (cov > max_a) {
	repeatP.resize(2*max_a+1);
	for(int a = max_a+1; a <= 2*max_a; a++) {
	  repeatP[a] = NULL;
	}
	max_a = 2*max_a;
      }
      if (max_cov < cov)
	max_cov = cov;
      repeatP[cov] = row;
    }
  }
  file.close();

  for(int a = 10; a <= max_cov; a++) {
    if (repeatP[a] == NULL)
      continue;
    for(int ell = a/2+1; ell < a; ell++) {
      double d_ell = (double) ell;
      double d_a = (double) a;
      double totP = 2.0*repeatP[a][1]*exp(-(d_ell-d_a/2.0)*(d_ell-d_a/2.0)*2.0/d_a);
      for(int k = 2; k < 5; k++) {
	if (ell < k/(k+1)*a) {
	  totP += repeatP[a][k];
	} else {
	  double d_k = (double)k;
	  totP += (d_k+1.0)*repeatP[a][k]*exp(-(d_k+1.0)/(2.0*d_a)*(d_ell-d_a*d_k/(d_k+1.0))*(d_ell-d_a*d_k/(d_k+1.0)));
	}
      }
      if (totP < epsilon) {
	th[a] = ell;
	prob[a] = totP;
	break;
      }
    }
  }

  for(int a = 0; a <= max_cov; a++) {
    if (repeatP[a] != NULL) {
      delete [] repeatP[a];
    }
  }
  
}
  

void Thresholds::writeThresholds(std::string filename)
{
  std::ofstream file(filename);

  for(int a = 1; a <= max_a; a++) {
    file << a << " " << th[a] << " " << prob[a] << "\n";
  }
  file.close();
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
      double p = std::stof(cell);
      while (a > max_a) {
	th.resize(max_a*2+1);
	for(int i = max_a+1; i <= max_a*2; i++) {
	  th[i] = i+1;
	}
	prob.resize(max_a*2+1);
	for(int i = max_a+1; i <= max_a*2; i++) {
	  prob[i] = 0.0;
	}
	max_a = 2*max_a;
      }
      th[a] = t;
      prob[a] = p;
    }
  }
  file.close();
}

