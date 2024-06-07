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

#ifndef PGMINFERENCE_H
#define PGMINFERENCE_H

#include <vector>
#include <map>
#include <list>
#include <set>

#include "factor.h"

// ============================================================================
// PGM INFERENCE
// ============================================================================

class PGMInference {
        typedef std::map<int, std::set<int> > AdjMatrix;
        typedef std::list<Factor> FactorList;

private:
        /**
         * Find a variable to eliminate based on the adjacency matrix. Selects
         * the variable (!= varToKeep) with the lowest number of neighbors
         * @param adjMatrix Adjacency matrix
         * @param varToKeep Variable that cannot be eliminated
         * @return Selected variable
         */
        static int findVarToEliminate(const AdjMatrix& adjMatrix, int varToKeep);

        /**
         * Given a list of factors, create a variable elimination ordering. Also
         * compute the number of variables in the largest intermediate factor
         * that would be produced when performing VE using this ordering.
         * @param factorList List of factors
         * @param varToKeep Variable that cannot be eliminated
         * @param elimOrder Vector with variables in elimination order
         * @return Number of values in largest intermediate factor
         */
        static size_t findElimOrdering(const FactorList& factorList,
                                       int varToKeep,
                                       std::vector<int>& elimOrder);

        /**
         * Eliminate a variable from a given a list of factors
         * @param factorList Input/output list of factors
         * @param varToEliminate Variable to eliminate
         */
        static void eliminateVar(FactorList& factorList,
                                 int varToEliminate);

public:
        /**
         * Given a list of factors as input, eliminate all variables with the
         * exception of the target variable. When successful the factorList
         * contains a single factor over the target variable
         * @param factorList List of input factors / output target factor
         * @param targetVar Target variable identifier
         * @param maxFactorSize Maximum number of values in intermediate factor
         * @return True VE completed within maxFactorSize constraints
         */
        static bool solveVE(FactorList& factorList, int targetVar,
                            size_t maxFactorSize);
};

#endif
