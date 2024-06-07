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

#ifndef FACTOR_H
#define FACTOR_H

#include <vector>
#include <cmath>
#include <iostream>

// ============================================================================
// ASSIGNMENT
// ============================================================================

class Assignment : public std::vector<int>
{
private:
        std::vector<int> card;          // cardinalities of variables
public:
        /**
         * Default constructor
         * @param card cardinalities of the variables
         */
        Assignment(const std::vector<int>& card) :
                std::vector<int>(card.size(), 0), card(card) {};

        /**
         * Increment the assignment in amortized constant time
         */
        void increment();
};

// ============================================================================
// INDEX CONVERTER
// ============================================================================

// Given two factors F_reference and F_subset where all variables in F_subset
// are also variables of the bigger factor F_reference. The index converter
// returns the index in F_subset while iterating through the index of
// F_reference. The conversion is done in (amortized) constant time.
class IndexConverter {

private:
        std::vector<int> refVar;        // variables in F_reference
        std::vector<int> refCard;       // cardinalities of F_reference
        std::vector<int> subVar;        // variables in F_subset
        std::vector<int> subCard;       // cardinalities of F_subset
        size_t refIndex, subIndex;      // index in F_reference and F_subset
        std::vector<int> refAss;        // assignment in F_reference
        std::vector<int> deltaSubIndex; // deltas to update F_subset

public:
        /**
         * Constructor
         * @param refVar List of reference variables
         * @param refCard List of reference cardinalities
         * @param subVar Subset of reference variables
         */
        IndexConverter(const std::vector<int>& refVar,
                       const std::vector<int>& refCard,
                       const std::vector<int>& subVar);

        /**
         * Get the cardinalities in F_subset
         * @return A const-ref to the cardinalities in F_subset
         */
        const std::vector<int>& getSubCard() const {
                return subCard;
        }

        /**
         * Get the current index of F_subset
         * @return Index of F_subset
         */
        size_t getSubIndex() const {
                return subIndex;
        }

        /**
         * Get the index of F_reference
         * @return Index of F_reference
         */
        size_t getRefIndex() const {
                return refIndex;
        }

        /**
         * Increment the index. In F_reference this simply moves to the next
         * index. In F_subset, the index may jump around.
         */
        void incrementIndex();
};

// ============================================================================
// FACTOR CLASS (with full table representation)
// ============================================================================

class Factor {
private:
        std::vector<int> var;           // variable IDs in the factor
        std::vector<int> card;          // cardinality of each variable
        std::vector<double> val;        // value of the assignments

        /**
         * Compute the sum of two numbers in log-space (numerically stable)
         * @param x First number (in log-space)
         * @param y Second number (in log-space)
         * @return Sum x+y (in log-space)
         */
        static double log_sum_exp(double x, double y) {
                if (x < y)
                        std::swap<double>(x, y);

                // for small values of exp(d), log(1+x) is almost equal to x
                double d = y - x;
                if (d < -10)    // looses no digits of accuracy
                        return x + exp(d);

                return x + log1p(exp(d));
        }

public:
        /**
         * Default constructor
         */
        Factor() {};

        /**
         * Create a full instantiated factor
         * @param var Vector with variable IDs (sorted)
         * @param card Vector with corresponding variable cardinalities
         * @param val Full table factor value
         */
        Factor(const std::vector<int>& var,
               const std::vector<int>& card,
               const std::vector<double>& val);

        /**
         * Create a full instantiated factor
         * @param var Vector with variable IDs (sorted)
         * @param card Vector with corresponding variable cardinalities
         * @param val Full table factor value
         */
        Factor(std::vector<int>&& var,
               std::vector<int>&& card,
               std::vector<double>&& val);

        /**
         * Get the number of variables in the factor
         * @return The number of variables in the factor
         */
        size_t getNumVar() const {
                return var.size();
        }

        /**
         * Get the variables in the factor
         * @return A const-reference to the variables in the factor
         */
        const std::vector<int>& getVar() const {
                return var;
        }

        /**
         * Set the variables in the factor
         * @param var Target variables
         */
        void setVar(const std::vector<int>& var) {
                this->var = var;
        }

        /**
         * Get the cardinalities in the factor
         * @return A const-reference to the cardinalities in the factor
         */
        const std::vector<int>& getCard() const {
                return card;
        }

        /**
         * Get the values in the factor
         * @return A const-reference to the values in the factor
         */
        const std::vector<double>& getVal() const {
                return val;
        }

        /**
         * Get the number of values in the factor
         * @return The number of values in the factor
         */
        size_t getNumVal() const {
                return val.size();
        }

        /**
         * Operator[] overloading
         * @return The value corresponding to the index
         */
        double operator[](size_t index) const {
                return val.at(index);
        }

        /**
         * Operator[] overloading
         * @return Reference to the value corresponding to the index
         */
        double& operator[](size_t index) {
                return val[index];
        }

        /**
         * Operator[] overloading
         * @return The value corresponding to the index
         */
        double operator[](const std::vector<int>& assignment) const {
                return val.at(assignmentToIndex(assignment));
        }

        /**
         * Operator[] overloading
         * @return Reference to the value corresponding to the index
         */
        double& operator[](const std::vector<int>& assignment) {
                return val[assignmentToIndex(assignment)];
        }

        /**
         * Convert an assignment to an index
         * @param assignment Factor assignment
         * @return Index
         */
        size_t assignmentToIndex(const std::vector<int>& assignment) const;

        /**
         * Convert an index to an assignment
         * @param index Index
         * @return Factor assignment
         */
        std::vector<int> indexToAssignment(size_t index) const;

        /**
         * Factor multiplication
         * @param rhs Right hand side value
         * @return A new factor containing the multiplied factors
         */
        Factor operator*(const Factor& rhs) const;

        /**
         * Marginalize a factor
         * @param varToKeep Variables to keep, other variables are marginalized
         * @return Marginal factor
         */
        Factor marginalize(const std::vector<int>& varToKeep) const;

        /**
         * Operator << overloading
         * @param os Output stream
         * @param F Factor F
         * @return Output stream with F appended to it
         */
        friend std::ostream& operator<<(std::ostream& os, const Factor& F);
        
         /**
         * Write out a factor block to a given file stream, according to the libDAI Factor graph format convention
         * Right now this function assumes there are no zero entries in the factor
         * @param file a .fg file to write to (has to be open)
         */
        void writeLibDAIFactorBlock(std::ofstream& file) const;
};

#endif
