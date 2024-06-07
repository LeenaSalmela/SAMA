#ifndef THRESHOLD_H
#define THRESHOLD_H

#include <vector>
#include <string>

#include "global.h"

class Thresholds
{
 private:
    std::vector<int> th;
    size_t max_a;
    
 public:
    Thresholds();

    int get(int a);

    void readThresholds(std::string filename);
};

#endif
