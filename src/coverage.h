/******************************************************************************
 *   Copyright (C) 2014 - 2022 Jan Fostier (jan.fostier@ugent.be)             *
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

#ifndef COVERAGE_H
#define COVERAGE_H

#include <cmath>
#include <iostream>
#include "util.h"

// ============================================================================
// CLASS PROTOTYPE
// ============================================================================

class Factor;

// ============================================================================
// MULTIPLICITY CLASS (Soft assignment)
// ============================================================================

class Multiplicity {

private:
        int firstMult;          // multiplicity correponding to P[0]
        std::vector<double> P;  // probabilities of multiplicities (log-space)

public:
        /**
         * Default constructor
         */
        Multiplicity() : firstMult(0) {
                P.push_back(0.0);       // log(1) = 0
        }

        /**
         * Constructor with single value
         * @param multiplicity Target multiplicity
         */
        Multiplicity(int multiplicity) : firstMult(multiplicity) {
                P.push_back(0.0);       // log(1) = 0
        }

        /**
         * Constructor from factor
         * @param firstMult First multiplicity in the factor
         * @param P Probabilities (log-space)
         */
        Multiplicity(int firstMult, const std::vector<double>& P) :
                firstMult(firstMult), P(P) {}

        /**
         * Get the probability corresponding to a particular multiplicity
         * @param mult Target multiplicity
         * @return Probability between [0..1]
         */
        double operator[](int mult) const;

        /**
         * Get the expected multiplicity (highest probability)
         * @return Expected multiplicity
         */
        int getExpMult() const;

        /**
         * Decrement the multiplicity (e.g. after some repeat resolution)
         * @param decreaseBy Value with which multiplicity should be decreased
         */
        void decrementMult(int decreaseBy = 1) {
                firstMult -= decreaseBy;
        }

        /**
         * Get the log(odds ratio) of the multiplicity being equal to the
         * expected multiplicity and the multiplicity being different from
         * the expected multiplicity. This provides a measure of confidence
         * in how certain we are that the assignment is correct
         * @return log( P(M = expMult) / P(M != expMult) )
         */
        double getExpMultLogOR() const;

        /**
         * Get the log probability of the multiplicity being equal to the
         * expected multiplicity
         * @return log( P(M = expMult) )
         */
        double getExpMultLProb() const;

        /**
         * Implicit conversion to integral type
         * @return Expected multiplicity
         */
        operator int() const {
                return getExpMult();
        }

        /**
         * Normalize the probabilities
         */
        void normalize();

        /**
         * Operator << overloading
         * @param out Output stream
         * @param m Multiplicity to display
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream &out, const Multiplicity &m);
};

// ============================================================================
// COVERAGE MODEL
// ============================================================================

class CovModel {

private:
        double errorLambda;     // mean of the error term
        double errorODF;        // error overdispersion factor

        double lambda;          // mean of multiplicity = 1
        double ODF;             // overdispersion factor
        std::vector<double> logw;       // log(weight) of the components

public:

        /**
         * Constructor
         * @param errorLambda Mean of the error term
         * @param lambda Mean of multiplicity = 1
         * @param w Weight of the components
         * @param ODF Overdispersion factor
         */
        CovModel(double errorLambda, double errorODF,
                 double lambda, double ODF,
                 const std::vector<double>& w);

        /**
         * Constructor from file
         * @param filename Input filename
         */
        CovModel(const std::string& filename);

        /**
         * Set a new value for lambda
         * @param newLambda target value for lambda
         */
        void setLambda(double newLambda) {
                lambda = newLambda;
        }

         /**
         * Set a new value for lambda_err
         * @param newErr target value for lambda_err
         */
        void setErrLambda(double newErr) {
                errorLambda = newErr;
        }

        /**
         * Get lambda_err
         * @return errorLambda
         */
        double getErrLambda() const {
                return errorLambda;
        }

        /**
         * Get the error overdispersion factor
         * @return Error overdispersion factor
         */
        double getErrorODF() const {
                return errorODF;
        }

        /**
         * Get lambda
         * @return lambda
         */
        double getLambda() const {
                return lambda;
        }

        /**
         * Get the overdispersion factor
         * @return The overdispersion factor
         */
        double getODF() const {
                return ODF;
        }

        /**
         * Get the weight for a given multiplicity
         * @param mult Multiplicity
         * @return The weight for the multiplicity
         */
        double getWeight(int mult) const {
                return (mult < (int)logw.size()) ?
                        exp(logw[mult]) : exp(logw.back());
        }

        /**
         * Set the weights (uniform accross multiplicities)
         * @param w Desired weight
         */
        void setWeightUniform(double w) {
                for (auto& it : logw)
                        it = log(w);
        }
        
        void setZeroWeight(double w) {
                logw[0] = log(w);
        }

        /**
         * Get the number of components in the model
         * @return The number of components
         */
        size_t getK() const {
                return logw.size();
        }

        /**
         * Get the hard assignment given the observed coverage
         * @param obsCov Observed coverage
         * @return The most likely multiplicity
         */
        int getExpMult(double obsCov) const;

        /**
         * Get the soft assignment given the observed coverage
         * @param obsCov Observed coverage
         * @param numAlt Number of alternatives
         * @return Multiplicity (soft assignment)
         */
        Multiplicity getMultSoft(double obsCov, int numAlt = 2) const;

        /**
         * Get the logprob of observing obsCov given the multiplicity
         * @param obsCov Observed coverage
         * @param mult Multiplicity (0 = error term)
         * @return P(obsCov | multiplicity)
         */
        double getLogProb(double obsCov, int mult) const;

        /**
         * Operator << overloading
         * @param out Output stream
         * @param c Coverage model to display
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream &out, const CovModel &c);

        /**
         * Write out parameters of coverage model
         * @param covFilename name of file to write to
         */
        void write(const std::string& covFilename) const;

        /**
         * Read parameters of coverage model from disk
         * @param covFilename name of file to read from
         */
        void read(const std::string& covFilename);

        /**
         * Write histogram and model fit to GNUplot file
         * @param baseFilename filename
         * @param hist Coverage data
         */
        void writeGnuplot(const std::string& baseFilename,
                          const std::map<int, double>& hist) const;

        /**
         * Get the coverage below which an error is most likely
         * @param epsilon Accuracy of the estimate
         * @return Coverage cutoff value
         */
        double getCovCutOff(double OR, double epsilon = 1e-3) const
        {
                // use bisection algorithm to find good estimate
                double left = errorLambda;
                double right = lambda;

                do {
                        double mid = 0.5 * (left + right);
                        double Pe = getLogProb(mid, 0);
                        double Pt = getLogProb(mid, 1);

                        if (Pe > OR * Pt)
                                left = mid;
                        else
                                right = mid;
                } while (right - left > epsilon);

                return 0.5 * (right + left);
        }
};

#endif
