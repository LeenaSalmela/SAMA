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

#ifndef UTIL_H
#define UTIL_H


#include <string>
#include <chrono>
#include <ctime>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <mutex>
#include <cassert>

// ============================================================================
// DEFINITIONS
// ============================================================================

#define MAX_TIMERS 16

// ============================================================================
// WORK LOAD BALANCER
// ============================================================================

class WorkLoadBalancer
{
private:
        size_t idxBegin;                // begin index to process
        size_t idxEnd;                  // end index to process
        size_t chunkSize;               // target chunk size
        size_t idxNext;                 // next node offset

        std::mutex mutex;               // mutex
        std::string message;            // message to print during progress

public:
        /**
         * Get a chunk of nodes (thread-safe)
         * @param chunkBegin Begin of chunk to handle
         * @param chunkEnd End of chunk to handle
         * @return False if there is no work left
         */
        bool getChunk(size_t& chunkBegin, size_t &chunkEnd);

        /**
         * Default constructor
         * @param idxBegin Begin index to process
         * @param idxEnd End index to process
         * @param chunkSize Maximum size per chunk
         * @param message Message to print during progress
         */
        WorkLoadBalancer(size_t idxBegin, size_t idxEnd, size_t chunkSize,
                         const std::string& message = "Processing") :
                idxBegin(idxBegin), idxEnd(idxEnd), chunkSize(chunkSize),
                idxNext(idxBegin), message(message) {}
};

// ============================================================================
// PHRED SCORE CONVERTER
// ============================================================================

class PhredConv
{
private:
        std::array<float, 256> lut;     // lookup table = P(correct)

public:
        /**
         * Default constructor (= use k-mer counts)
         */
        PhredConv() {
                for (int i = 0; i < 256; i++)
                        lut[i] = 1.0f;
        }

        /**
         * Initialize lookup table with the Phred score
         * @param base ASCII value corresponding to Q = 0
         */
        void enablePhred(int base = 33) {
                for (int i = base; i < 256; i++)
                        lut[i] = 1.0 - pow(10.0, -(i-base)/10.0);
        }

        /**
         * Get P that nucleotide is correct given its ASCII Phred score
         * @param c Phred score encoded as ASCII character
         * @return Probability that nucleotide is correct
         */
        float getProb(char c) {
                return lut[(int)c];
        }

        /**
         * Get Phred score encoded as ASCII character
         * @param p Probability that nucleotide is incorrect
         * @return Phred score encoded as ASCII character
         */
        char getPhred(double p, int base=33) {
	  if (-10.0*log10(p) > 40) {
	    return base+40;
	  } else {
	    return (char)(-10.0*log10(p)+base);
	  }
        }


  
};

// ============================================================================
// UTILITY CLASS WITH DIVERSE AUXILIARY ROUTINES
// ============================================================================

class Util
{
private:
        static double prevProgress;
        static int currentTimer;
        static std::chrono::time_point<std::chrono::system_clock> startTime[MAX_TIMERS];
        static PhredConv phredConv;

public:
        /**
         * Create a string with a human readable version of a time period
         * @param time Time period (expressed in s)
         * @return String with a human readable version of a time period
         */
        static std::string humanTime(double time);

        /**
         * Create a string with a human readable version of a (file) size
         * @param size Size (expressed in bytes)
         * @return String with a human readable version of a (file) size
         */
        static std::string humanSize(size_t size);

        /**
         * Write a progress indicator to stdout
         * @param str String to print before percentage
         * @param curr Current progress thus far
         * @param max Maximum progress
         */
        static void progress(const std::string& str, double curr, double max);

        /**
         * Write a progress indicator to stdout followed by "(100.0%)"
         * @param str String to print before percentage
         * @param elapsed Elapsed amount of time (optional)
         */
        static void progressEnd(const std::string& str, double elapsed = -1.0);

        /**
         * Start a chronometer
         */
        static void startChrono();

        /**
         * Stop the chronometer
         * @return The time in
         */
        static double stopChrono();

        /**
         * Stop the chronometer and return a human readable string
         * @return A human readable string containg the elapsed time
         */
        static std::string stopChronoStr() {
                return humanTime (stopChrono());
        }

        /**
         * Get a string containing the date and time
         * @return string containing date and time
         */
        static std::string getDateTime();

        /**
         * Compute the specificity
         * @param TN True negatives
         * @param FP False positives
         * The specificity
         */
        static double getSpecificity(double TN, double FP) {
                return ((TN+FP) == 0) ? 1.0 : TN / (TN + FP);
        }

        /**
         * Compute the sensitivity
         * @param TP True positives
         * @param FN False negatives
         * The sensitivity
         */
        static double getSensitivity(double TP, double FN) {
                return ((TP+FN) == 0) ? 1.0 : TP / (TP + FN);
        }

        /**
         * Check wether a file exists
         * @param filename
         * @return True of false
         */
        static bool fileExists(const std::string& filename) {
                std::ifstream file(filename.c_str(), std::ios::in);
                bool OK = file.good();
                file.close();
                return OK;
        }

        /**
         * Compute the relative distance between two numbers
         * @param a number one
         * @param b number two
         * @return abs((a-b)/min(a,b)) OR 0 when both a=0 and b=0
         */
        static double relDiff(double a, double b) {
                if (a == 0.0 && b == 0.0)
                        return 0.0;
                return std::max(std::abs((a - b) / a), std::abs((a - b) / b));
        }

        /**
         * Compute the probability p(k) from a Poisson distribution with mean mu
         * @param k Number of observations
         * @param mu Average
         * @return The probability p(k)
         */
        static double poissonPDF(unsigned int k, double mu);

        /**
         * Compute the logprob p(k) from a Poisson distribution with mean mu
         * @param k Number of observations
         * @param mu Average
         * @return The probability p(k)
         */
       	static double logPoissonPDF(unsigned int k, double mu);
        static double logPoissonPDF(double k, double mu);

        /**
         * Compute the probability ratio p(k, mu1) / p(k, mu2)
         * @param k Number of observations
         * @param mu1 Expected number of observation (mean of distribution)
         * @param mu2 Expected number of observation (mean of distribution)
         * @return The probability p(k)
         */
        static double poissonPDFratio(unsigned int k, double mu1, double mu2);

        /**
         * Compute the probability p(k) from a negative bionomial(mu, var)
         * @param k Number of observations
         * @param mu Average
         * @param var Variance
         * @return The probability p(k)
         */
        static double negbinomialPDF(unsigned int k, double mu, double var);

        /**
         * Compute the log(probability p(k)) from a negative bionomial(mu, var)
         * @param k Number of observations
         * @param mu Average
         * @param var Variance
         * @return The probability p(k)
         */
        static double logNegbinomialPDF(unsigned int k, double mu, double var);
        static double logNegbinomialPDF(double k, double mu, double var);

        /**
         * Compute the probability ratio p(k, mu1) / p(k, mu2)
         * @param k Number of observations
         * @param mu1 Mean of the distribution
         * @param var1 Variance of the first distribution
         * @param mu2 Mean of the second distribution
         * @param var2 Variance of the second distribution
         * @return The probability p(k)
         */
        static double negbinomialPDFratio(unsigned int k,
                                          double mu1, double var1,
                                          double mu2, double var2);

        /**
         * Compute the probability p(k) for the geometric distribution
         * @param k Number of observations
         * @param mu Mean of the distribution
         * @param mu2 Variance of the distribution
         * @return The probability p(k)
         */
        static double geometricPDF(unsigned int k, double mu);

        /**
         * Compute the log(probability p(k)) for the geometric distribution
         * @param k Number of observations
         * @param mu Mean of the distribution
         * @param mu2 Variance of the distribution
         * @return The probability p(k)
         */
        static double logGeometricPDF(unsigned int k, double mu);

        /**
         * Compute the probability ratio p(mu1) / p(k, mu2)
         * @param k Number of observations
         * @param mu1 Mean of the geometric distribution
         * @param mu2 Mean of the negative bionomial distribution
         * @param var2 Variance of the negative bionomial distribution
         * @return The probability p(k)
         */
        static double geometricnegbinomialPDFratio(unsigned int k, double mu1,
                                                   double mu2, double var2);

        /**
         * Compute the percentage of two size_t numbers
         * @param nom Nominator
         * @param den Denominator
         * @return The percentage
         */
        static double toPercentage(size_t nom, size_t den) {
                if (den == 0)
                        return 0;
                return 100.0 * double(nom) / double(den);
        }

        /**
         * Fit a negative bionomial to truncated data
         * @param data Data to which the model is fitted < x, y >
         * @param mu Mean of NB
         * @param ODF Overdispersion factor of NB (= variance / mu)
         * @param w Weight of NB
         * @param epsilon Convergence criterion
         * @param maxIter Maximum number of iterations
         * @return Number of iterations (maxIter + 1 if not converged)
         */
        static int fitTruncNegBinomEM(std::map<unsigned int, double>& data,
                                      double& mu, double& ODF,
                                      double& w, double epsilon, int maxIter);

        /**
         * Fit normal mixture model to data points
         * @param data Data to which the model is fitted < x, y >
         * @param mu Means of model components
         * @param sigma2 Variance of model components
         * @param MC Weights of models components
         * @param maxIteration Maximum number of iterations
         */
        static void binomialMixtureEM(const std::map<unsigned int, double>& data,
                                      std::vector<double>& mu,
                                      std::vector<double>& sigma2,
                                      std::vector<double>& MC,
                                      int maxIterations = 20);

        /**
         * Compute the sum of two numbers in log-space (numerically stable)
         * @param x First number (in log-space)
         * @param y Second number (in log-space)
         * @return Sum x+y (in log-space)
         */
        static double log_sum_exp(double x, double y);

        static void enablePhred(int base = 33) {
                phredConv.enablePhred(base);
        }

        static void writeSeqWrap(std::ostream& ofs, const std::string& s,
                                 size_t wrap) {
                for (size_t c = 0; c < s.size(); c += wrap)
                        ofs << s.substr(c, wrap) << "\n";
        }

        static void writeProbWrap(std::ostream& ofs, const std::vector<double>& p, size_t wrap) {
                char *buf = new char[wrap+1];
	  
                for (size_t c = 0; c < p.size(); c +=wrap) {
		  for(size_t i = 0; i < wrap && c+i < p.size(); i++) {
		    buf[i] = phredConv.getPhred(p[c+i]);
		  }
		  if (c+wrap <= p.size())
		    buf[wrap] = '\0';
		  else
		    buf[p.size() % wrap] = '\0';
		  ofs << buf << "\n";
		}
        }

        /**
         * Compute P that sequence is correct given ASCII phred scores
         * @param phred Phred scores encoded as ASCII string
         * @param b Begin position
         * @param e End position
         * @return P([b, e[ is correct) = product of individual probabilities
         */
        static double phred2prob(const std::string& phred, size_t b, size_t e);

        /**
         * Compute the factorial of a number
         * @param n Number (>= 0)
         * @return The factorial product
         */
        static size_t factorial(int n) {
                assert(n >= 0);
                size_t res = 1;
                for (int i = 2; i <= n; i++)
                        res *= i;
                return res;
        }

        /**
         * Compute the ratio a! / b!
         * @param a Number a (>= 0)
         * @param b Number b (>= 0) AND (a >= b OR a = 0; b = 1)
         * @return The ratio a! / b!
         */
        static size_t factRatio(int a, int b) {
                assert((a >= 0) && (b >= 0));
                assert((a == 0 && b == 1) || (a >= b));
                size_t res = 1;
                for (int i = b+1; i <= a; i++)
                        res *= i;
                return res;
        }
};

#endif
