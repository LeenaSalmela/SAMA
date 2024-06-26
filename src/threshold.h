#ifndef THRESHOLD_H
#define THRESHOLD_H

#include <vector>
#include <string>

#include "global.h"

class Thresholds
{
 private:
    std::vector<int> th;
    std::vector<double> prob;
    size_t max_a;
    
 public:
    Thresholds();

    int get(int a);

  void computeThresholds(std::string filename, double epsilon);
  
    void writeThresholds(std::string filename);

    void readThresholds(std::string filename);
};

#endif
