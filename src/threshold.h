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
    std::vector<double *> totP;
    size_t max_a;
    
 public:
    Thresholds();

    int get(int a);
    double getProb(int nodeCount, int arcCount);
  
    void computeThresholds(std::string filename, double epsilon);
  
    void writeThresholds(std::string thresholdfile, std::string probfile);

    void readThresholds(std::string thresholdfile, std::string probfile);
};

#endif
