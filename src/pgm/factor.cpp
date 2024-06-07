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

#include <cassert>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cmath>
#include <fstream>

#include "factor.h"

using namespace std;

// ============================================================================
// AUXILIARY ROUTINES
// ============================================================================

/**
 * Returns indexes [0 1 2 ...] sorted according to values in v
 * @param v Vector of values
 * @return Sorted indices
 */
template <typename T>
vector<size_t> sort_indexes(const vector<T>& v)
{
        // initialize original index locations
        vector<size_t> idx(v.size());
        iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        sort(idx.begin(), idx.end(),
             [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

        return idx;
}

// ============================================================================
// ASSIGNMENT
// ============================================================================

void Assignment::increment()
{
        // amortized constant-time procedure to update subIndex:
        // find the most significant variable position j in the reference
        // assignment that is updated. This means all variables at positions
        // smaller than j are reset to zero.
        for (size_t j = 0; j < size(); j++) {
                (*this)[j]++;
                if ((*this)[j] != card[j])
                        break;
                (*this)[j] = 0;
        }
}

// ============================================================================
// INDEX CONVERTOR
// ============================================================================

IndexConverter::IndexConverter(const vector<int>& refVar,
                               const vector<int>& refCard,
                               const vector<int>& subVar) :
                               refVar(refVar), refCard(refCard), subVar(subVar),
                               refIndex(0), subIndex(0),
                               refAss(refVar.size(), 0)
{
        // sort the reference variables
        vector<int> refVarSort = refVar;
        sort(refVarSort.begin(), refVarSort.end());
        vector<size_t> refIdxSrt2Unsrt = sort_indexes(refVar);

        // make sure there are no doubles in the variable names
        assert(adjacent_find(refVarSort.begin(), refVarSort.end()) == refVarSort.end());

        // sort the subfactor variables
        vector<int> subVarSort = subVar;
        sort(subVarSort.begin(), subVarSort.end());
        vector<size_t> subIdxSrt2Unsrt = sort_indexes(subVar);

        // make sure there are no doubles in the variable names
        assert(adjacent_find(subVarSort.begin(), subVarSort.end()) == subVarSort.end());

        // make sure all variables in subVar are contained within refVar
        assert(includes(refVarSort.begin(), refVarSort.end(),
                        subVarSort.begin(), subVarSort.end()));

        // ref2sub translates an index in refVar to the corresponding index in subVar
        vector<int> ref2sub(refVar.size(), -1);
        // also populate the subCard variable
        subCard.resize(subVar.size());
        for (size_t is = 0, ir = 0; is < subVar.size(); is++) {
                while (refVarSort[ir] != subVarSort[is])
                        ir++;
                assert(refVar[refIdxSrt2Unsrt[ir]] == subVar[subIdxSrt2Unsrt[is]]);
                ref2sub[refIdxSrt2Unsrt[ir]] = subIdxSrt2Unsrt[is];
                subCard[subIdxSrt2Unsrt[is]] = refCard[refIdxSrt2Unsrt[ir]];
        }

        // compute the multipliers in F_subset
        vector<int> subMult(subVar.size(), 1);
        for (size_t i = 1; i < subVar.size(); i++)
                subMult[i] = subMult[i-1] * subCard[i-1];

        // deltaSubIndex[i] contains the values that need to be added to
        // subIndex when variable at position i the most significant value
        // in the reference assignment that is updated during incrementIndex
        deltaSubIndex.resize(refVar.size(), 0);
        int reset_ctr = 0;
        for (size_t i = 0; i < refVar.size(); i++) {
                // incrementing the variable at position i will reset all
                // variables 0..i-1 back to zero, causing a decrease of subIndex
                deltaSubIndex[i] -= reset_ctr;

                // if the variable i is not in subVar, we're done
                if (ref2sub[i] == -1)
                        continue;

                // incrementing variable position i will increase subIndex
                deltaSubIndex[i] += subMult[ref2sub[i]];
                reset_ctr += (refCard[i]-1) * subMult[ref2sub[i]];
        }
}

void IndexConverter::incrementIndex()
{
        refIndex++;     // simply move to the next index in refIndex

        // amortized constant-time procedure to update subIndex:
        // find the most significant variable position j in the reference
        // assignment that is updated. This means all variables at positions
        // smaller than j are reset to zero.
        for (size_t j = 0; j < refAss.size(); j++) {
                refAss[j]++;
                if (refAss[j] == refCard[j]) {
                        refAss[j] = 0;
                } else {
                        // we're at the most significant variable position j
                        subIndex += deltaSubIndex[j];
                        break;
                }
        }
}

// ============================================================================
// FACTOR CLASS (with full table representation)
// ============================================================================

Factor::Factor(const vector<int>& var,
               const vector<int>& card,
               const vector<double>& val) :
               var(var), card(card), val(val)
{
        // make sure var and card have equal size
        assert(var.size() == card.size());

        // make sure all values of card are > 0
        assert(find_if(card.begin(), card.end(), [](int i){return i < 1;}) == card.end());

#ifdef DEBUG
        // make sure there are no doubles in var
        vector<int> temp = var;
        sort(temp.begin(), temp.end());
        assert(adjacent_find(temp.begin(), temp.end()) == temp.end());

        // make sure val has the correct size (product of cardinalities)
        size_t numVal = accumulate(card.begin(), card.end(), 1ull, multiplies<size_t>());
        assert(val.size() == numVal);
#endif
}

Factor::Factor(vector<int>&& var,
               vector<int>&& card,
               vector<double>&& val) :
               var(var), card(card), val(val)
{
        // make sure var and card have equal size
        assert(var.size() == card.size());

        // make sure all values of card are > 0
        assert(find_if(card.begin(), card.end(), [](int i){return i < 1;}) == card.end());

#ifdef DEBUG
        // make sure there are no doubles in var
        vector<int> temp = var;
        sort(temp.begin(), temp.end());
        assert(adjacent_find(temp.begin(), temp.end()) == temp.end());

        // make sure val has the correct size (product of cardinalities)
        size_t numVal = accumulate(card.begin(), card.end(), 1ull, multiplies<size_t>());
        assert(val.size() == numVal);
#endif
}

size_t Factor::assignmentToIndex(const vector<int>& assignment) const
{
        assert(assignment.size() == var.size());

        size_t index = 0;
        for (size_t i = 0, multiplier = 1; i < assignment.size(); i++) {
                index += assignment[i] * multiplier;
                multiplier *= card[i];
        }

        assert(index < val.size());
        return index;
}

vector<int> Factor::indexToAssignment(size_t index) const
{
        assert(index < val.size());

        vector<int> multiplier(card.size(), 1);
        for (size_t i = 1; i < card.size(); i++)
                multiplier[i] = multiplier[i-1] * card[i-1];

        vector<int> assignment(card.size());
        for (size_t i = assignment.size(); i-- > 0; ) {
                assignment[i] = index / multiplier[i];
                index -= assignment[i] * multiplier[i];
        }

        return assignment;
}

Factor Factor::operator*(const Factor& rhs) const
{
        // treat empty factors as identity factors
        if (rhs.var.empty())
                return *this;
        if (var.empty())
                return rhs;

        // create a vector with joint variables and cardinalities
        vector<pair<int, int> > merge;
        merge.reserve(var.size() + rhs.var.size());
        for (size_t i = 0; i < var.size(); i++)
                merge.push_back(make_pair(var[i], card[i]));
        for (size_t i = 0; i < rhs.var.size(); i++)
                merge.push_back(make_pair(rhs.var[i], rhs.card[i]));

        // sort the vector and remove all doubles
        sort(merge.begin(), merge.end());
        auto last = unique(merge.begin(), merge.end());
        merge.erase(last, merge.end());

        // create the merged variable and cardinalities
        vector<int> newVar, newCard;
        newVar.reserve(merge.size());
        newCard.reserve(merge.size());
        for (auto it : merge) {
                newVar.push_back(it.first);
                newCard.push_back(it.second);
        }

        // assert there are no double in varM. This guarantees that shared
        // variables between *this and rhs have the same cardinality
        assert(unique(newVar.begin(), newVar.end()) == newVar.end());

        // idxConv will convert an index in the original factor to a
        // corresponding index in the marginal factor
        IndexConverter idxConv1(newVar, newCard, var);
        IndexConverter idxConv2(newVar, newCard, rhs.var);

        // compute the factor product
        int numVal = accumulate(newCard.begin(), newCard.end(), 1, multiplies<int>());
        vector<double> newVal(numVal);

        for (size_t i = 0; i < newVal.size(); i++) {
                size_t idx1 = idxConv1.getSubIndex();
                size_t idx2 = idxConv2.getSubIndex();
                newVal[i] = val[idx1] + rhs.val[idx2];

                idxConv1.incrementIndex();
                idxConv2.incrementIndex();
        }

        return Factor(move(newVar), move(newCard), move(newVal));
}

Factor Factor::marginalize(const vector<int>& margVar) const
{
        // idxConv will convert an index in the original factor to a
        // corresponding index in the marginal factor
        IndexConverter idxConv(var, card, margVar);

        // get the new cardinalities
        vector<int> margCard = idxConv.getSubCard();

        // get the new values
        int numMargVal = accumulate(margCard.begin(), margCard.end(),
                                    1, multiplies<int>());
        // initialize to log(0) = -inf -> initialize to smallest possible double
        vector<double> margVal(numMargVal, numeric_limits<double>::lowest());

        for (size_t i = 0; i < val.size(); i++) {
                size_t j = idxConv.getSubIndex();
                margVal[j] = log_sum_exp(margVal[j], val[i]);
                idxConv.incrementIndex();
        }

        vector<int> copyVar(margVar);
        return Factor(move(copyVar), move(margCard), move(margVal));
}

ostream& operator<<(ostream& os, const Factor& F) {
        os << "Variables in F: ";
        for (size_t i = 0; i < F.var.size(); i++)
                os << F.var[i] << "\t";
        os << "\nCardinalities of variables: ";
        for (size_t i = 0; i < F.card.size(); i++)
                os << F.card[i] << "\t";
        os << "\nValues of F:\n";
        Assignment assignment(F.card);
        for (size_t i = 0; i < F.val.size(); i++) {
                for (size_t j = 0; j < assignment.size(); j++)
                        os << assignment[j] << "\t";
                os << ": " << F.val[i] << "\n";
                assignment.increment();
        }

        return os;
}

void Factor::writeLibDAIFactorBlock(std::ofstream& file) const{
        if( ! file.is_open())
                return;
        
        file << "\n" << getNumVar() << "\n";
        for (size_t i = 0; i < var.size(); i++)
                file << var[i] << " ";
        file << "\n";
        for (size_t i = 0; i < card.size(); i++)
                file << card[i] << " ";
        file << "\n" << val.size() << "\n";
        for (size_t i = 0; i < val.size(); i++) 
            file << i << " " << exp(val[i]) << "\n";
}
