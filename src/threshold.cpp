/******************************************************************************
 *   Copyright (C) 2024 Leena Salmela (leena.salmela@helsinki.fi)             *
 *   This file is part of SAMA                                                *
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Affero General Public License as published *
 *   by the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU Affero General Public License *
 *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/
#include "threshold.h"

#include <iostream>
#include <fstream>
#include <tuple>
#include <math.h>
#include <stdlib.h>
#include <string.h>

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
  totP.resize(max_a+1);
  for(int i = 0; i <= max_a; i++) {
    totP[i] = NULL;
  }
}

int Thresholds::get(int a)
{
  if (a > max_a) {
    return a+10000;
  }

  return th[a];
}

double Thresholds::getProb(int nodeCount, int arcCount) {

  if (nodeCount > max_a) {
    return 1.0;
  }

  if (arcCount-nodeCount/2-1 > nodeCount/2) {
    return 1.0;
  }

  if (arcCount-nodeCount/2-1 < 0) {
    return 0.0;
  }

  return totP[nodeCount][arcCount-nodeCount/2-1];
}
  

// Compute total probabilities for each [node_coverage, edge_coverage] pair
// Deduce thresholds
void Thresholds::computeThresholds(std::string filename, double epsilon) {
  std::ifstream file(filename);
  std::string cell;

  std::vector<double *> repeatP;
  int max_cov = 10;
  int est_cov = -1.0;
  int min_cov = -1;
  double est_cov_maxP = -1.0;
  
  std::cout << epsilon << std::endl;
  
  repeatP.resize(max_a+1);
  for(int a = 0; a <= max_a; a++) {
    repeatP[a] = NULL;
  }

  bool min_cov_set = false;
  while(!file.eof())
  {
    if (std::getline(file, cell, '\t'))
    {	
      double *row = new double[5];
      int cov = std::stoi(cell);
      std::getline(file, cell, '\t');  // 0 copies
      double p0 = std::strtod(cell.c_str(), NULL);
      double s = 0.0;
      for(int r = 1; r <= 5; r++) {
	std::getline(file, cell, '\t');  // r copies
	double p = std::strtod(cell.c_str(), NULL);
	row[r-1] = p;
	s+= p;
      }
      std::getline(file, cell, '\n');
      if (est_cov_maxP < row[0]) {
	est_cov = cov;
	est_cov_maxP = row[0];
      }
      if (p0 > row[0] && ! min_cov_set)
	min_cov = cov+1;
      else
	min_cov_set = true;
      for(int r = 1; r <= 5; r++) {
	row[r-1] = row[r-1]/s;
      }
      while (cov > max_a) {
	repeatP.resize(2*max_a+1);
	th.resize(2*max_a+1);
	prob.resize(2*max_a+1);
	totP.resize(2*max_a+1);
	for(int a = max_a+1; a <= 2*max_a; a++) {
	  repeatP[a] = NULL;
	  th[a] = a+10000;
	  prob[a] = 0.0;
	  totP[a] = NULL;
	}
	max_a = 2*max_a;
      }
      if (max_cov < cov)
	max_cov = cov;
      repeatP[cov] = row;
    }
  }
  file.close();

  if (min_cov < 0)
    min_cov = 5;
  
  std::cout << "Minimum coverage: " << min_cov << std::endl;
  std::cout << "Estimated coverage: " << est_cov << std::endl;
  
  
  for(int a = 1; a <= max_cov; a++) {
    bool threshold_found = false;
    if (repeatP[a] == NULL)
      continue;
    totP[a] = new double[a/2+1];
    //std::cout << "a: " << a << std::endl;
    for(int ell = a/2+1; ell <= a; ell++) {
      //std::cout << "  ell: " << ell << std::endl;
      double d_ell = (double) ell;
      double d_a = (double) a;
      //totP[a][ell-a/2-1] = 2.0*repeatP[a][1]*exp(-(d_ell-d_a/2.0)*(d_ell-d_a/2.0)*2.0/d_a);  // Chernoff
      totP[a][ell-a/2-1] = 0.0;
      
      //std::cout << "    totP: " << totP[a][ell-a/2-1] << std::endl;
      for(int k = 1; k < 5; k++) {
	double d_k = (double)k;
	if (ell <= d_k/(d_k+1.0)*a) {
	  totP[a][ell-a/2-1] += repeatP[a][k];
	} else {
	  //totP[a][ell-a/2-1] += (d_k+1.0)*repeatP[a][k]*exp(-(d_k+1.0)/(2.0*d_a)*(d_ell-d_a*d_k/(d_k+1.0))*(d_ell-d_a*d_k/(d_k+1.0))); // Chernoff
	  
	  if (a == ell) {
	    totP[a][ell-a/2-1] += (d_k+1.0)*repeatP[a][k]*pow((d_k/(d_k+1)), a);
	  } else {
	    double aa = d_ell / d_a;
	    double p = d_k/(d_k+1.0);
	    double r = p*(1.0-aa) / (aa*(1.0-p));
	    double Dap = aa*log(aa/p) + (1.0-aa)*log ((1.0-aa)/(1.0-p));
	    totP[a][ell-a/2-1] += (d_k+1.0)*repeatP[a][k]*1.0/(1.0-r)* 1.0 / sqrt(2.0*M_PI*aa*(1.0-aa)*d_a) * exp(-d_a*Dap);
	  }
	}
	//std::cout << "    totP: " << totP[a][ell-a/2-1] << std::endl;
      }
      // totP is an estimate so it can be >1, reset to 1 if so
      if (totP[a][ell-a/2-1] > 1.0) {
	totP[a][ell-a/2-1] = 1.0;
      }
      if (!threshold_found && totP[a][ell-a/2-1] < epsilon) {
	th[a] = ell;
	prob[a] = totP[a][ell-a/2-1];
	//std::cout << "    Threshold!" << std::endl;
	threshold_found = true;
      }
    }
  }

  for(int a = 0; a <= max_cov; a++) {
    if (repeatP[a] != NULL) {
      delete [] repeatP[a];
    }
  }
  
}
  

void Thresholds::writeThresholds(std::string thresholdfile, std::string probfile)
{
  std::ofstream file(thresholdfile);

  for(int a = 1; a <= max_a; a++) {
    file << a << " " << th[a] << " " << prob[a] << "\n";
  }
  file.close();

  std::ofstream pfile(probfile);

  for(int a = 1; a <= max_a; a++) {
    if (totP[a] != NULL) {
      pfile << a;
      for(int ell = a/2+1; ell <= a; ell++) {
	pfile << " " << totP[a][ell-a/2-1];
      }
      pfile << "\n";
    }
  }
  pfile.close();

}


void Thresholds::readThresholds(std::string thresholdfile, std::string probfile)
{
  std::ifstream file(thresholdfile);
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
	totP.resize(max_a*2+1, NULL);
	for(int i = max_a+1; i <= max_a*2; i++) {
	  totP[i] = NULL;
	}
	max_a = 2*max_a;
      }
      th[a] = t;
      prob[a] = p;
    }
  }
  file.close();

  std::ifstream pfile(probfile);
  
  while(!pfile.eof())
  {
    if (std::getline(pfile, cell, ' '))
    {	
      int a = std::stoi(cell);

      if (a > max_a) {
	std::cout << "Incompatible threshold and probability files" << std::endl;
	exit(2);
      }
      totP[a] = new double[a/2+1];
      std::getline(pfile, cell);
      char *buf = new char[cell.length()+1];
      strncpy(buf, cell.c_str(), cell.length()+1);
      buf[cell.length()] = '\0';
      char *t = strtok(buf, " ");
      int ell = 0;
      while(t!= NULL) {
	totP[a][ell] = atof(t);
	ell++;
	t = strtok(NULL, " ");
      }
      delete [] buf;
    }
  }
  pfile.close();


}



