/******************************************************************************
 *   Copyright (C) 2018 - 2022 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of Detox                                               *
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

#include <algorithm>
#include <cassert>

#include "pgminference.h"

using namespace std;

// ============================================================================
// PGM INFERENCE
// ============================================================================

int PGMInference::findVarToEliminate(const AdjMatrix& adjMatrix, int varToKeep)
{
        // find the random variable that has the fewest neighbors
        unsigned int minNeighbors = adjMatrix.size();
        int varToEliminate = 0;

        for (const auto& it : adjMatrix) {
                const int& var = it.first;
                const set<int>& neighbors = it.second;

                if (var == varToKeep)
                        continue;

                if (neighbors.size() < minNeighbors) {
                        minNeighbors = neighbors.size();
                        varToEliminate = var;
                }
        }

        return varToEliminate;
}

void PGMInference::eliminateVar(FactorList& factorList, int varToEliminate)
{
        // pointer to the first factor that has varToEliminate in its scope
        Factor* prod = NULL;

        // multiply all factors that contain variable varToEliminate
        for (auto it = factorList.begin(); it != factorList.end(); ) {
                Factor& F = *it;
                const vector<int>& var = F.getVar();
                if (find(var.begin(), var.end(), varToEliminate) == var.end()) {
                        it++;
                } else {        // F contains varToEliminate
                        if (prod == NULL) {     // first encounter?
                                prod = &F;
                                it++;
                        } else {
                                *prod = *prod * F;
                                it = factorList.erase(it);
                        }
                }
        }

        // find out what variables are contained in the product
        vector<int> varToKeep = prod->getVar();
        varToKeep.erase(remove(varToKeep.begin(), varToKeep.end(),
                               varToEliminate), varToKeep.end());
        *prod = prod->marginalize(varToKeep);
}

size_t PGMInference::findElimOrdering(const FactorList& factorList,
                                      int varToKeep, vector<int>& elimOrder)
{
        // create the adjacency matrix for the random variables
        AdjMatrix adjMatrix;
        for (const auto& F : factorList) {
                const vector<int>& varF = F.getVar();
                for (int var1 : varF)
                        for (int var2 : varF)
                                if (var1 != var2)
                                        adjMatrix[var1].insert(var2);
        }

        // handle the case where there is only one variable
        if (adjMatrix.empty())
                return 0;

        // get the cardinality of all the variables in the factorList
        map<int, int> var2card;
        for (const auto& F : factorList) {
                const vector<int>& varF = F.getVar();
                const vector<int>& cardF = F.getCard();

                for (size_t i = 0; i < varF.size(); i++)
                        var2card[varF[i]] = cardF[i];
        }

        // make sure the target variable is contained in the factorList
        assert(adjMatrix.find(varToKeep) != adjMatrix.end());

        // now find an elimination ordering
        elimOrder.clear();
        elimOrder.reserve(adjMatrix.size() - 1);

        size_t largestFactorSize = 1;
        while (adjMatrix.size() > 1) {  // eliminate all but one variables
                // find the best variable to eliminate
                int varToEliminate = findVarToEliminate(adjMatrix, varToKeep);
                elimOrder.push_back(varToEliminate);

                // find the size of the intermediate factor at this step
                size_t thisFactorSize = var2card[varToEliminate];
                const set<int>& nb = adjMatrix[varToEliminate];
                for (int var : nb)
                        thisFactorSize *= (size_t)var2card[var];
                largestFactorSize = max<size_t>(thisFactorSize, largestFactorSize);

                // update the connections in the adjacency matrix
                for (int var1 : nb) {
                        for (int var2 : nb)
                                if (var1 != var2)
                                        adjMatrix[var1].insert(var2);
                        // remove the connection to varToEliminate
                        adjMatrix[var1].erase(varToEliminate);
                }

                adjMatrix.erase(varToEliminate);
        }

        return largestFactorSize;
}

bool PGMInference::solveVE(FactorList& factorList, int targetVar,
                           size_t maxFactorSize)
{
        // get out when the factor list is empty
        if (factorList.empty())
                return false;

        // construct an elimination ordering and check feasibility
        vector<int> elimOrder;

        size_t factorSize = findElimOrdering(factorList, targetVar, elimOrder);
        if (factorSize > maxFactorSize)
                return false;

        // eliminate all variables
        for (int varToEliminate : elimOrder)
                eliminateVar(factorList, varToEliminate);

        // multiply remaining (singleton) factors
        while (factorList.size() > 1) {
                factorList.front() = factorList.front() * factorList.back();
                factorList.pop_back();
        }

        return true;
}
