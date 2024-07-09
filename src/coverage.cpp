/******************************************************************************
 *   Copyright (C) 2024 Leena Salmela (leena.salmela@helsinki.fi)             *
 *   This file has been modified for MAGA                                     *
 *                                                                            *
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

#include "dbgraph.h"
#include "kmernpp.h"
#include "settings.h"
#include "coverage.h"

using namespace std;

// ============================================================================
// MULTIPLICITY CLASS
// ============================================================================

int Multiplicity::getExpMult() const
{
        // sanity check
        assert(!P.empty());

        // return the multiplicity with the highest probability
        int maxOffset = 0; double maxProb = P[0];
        for (int i = 1; i < (int)P.size(); i++) {
                if (P[i] > maxProb) {
                        maxProb = P[i];
                        maxOffset = i;
                }
        }

        return firstMult + maxOffset;
}

double Multiplicity::getExpMultLogOR() const
{
        // if we only have one option, the odds ratio is infinite
        if (P.size() <= 1)
                return numeric_limits<double>::max();

        // get the most likely option
        int expMult = getExpMult();

        // compute the sum of the probabilities of the other options (log-space)
        double logSum = numeric_limits<double>::lowest();       // log(0) = -oo
        for (int i = 0; i < (int)P.size(); i++) {
                if (i + firstMult == expMult)
                        continue;
                logSum = Util::log_sum_exp(logSum, P[i]);
        }

        return P[expMult-firstMult] - logSum;
}

double Multiplicity::getExpMultLProb() const
{
        // if we only have one option, the log of the probability is zero
        if (P.size() <= 1)
                return 0.0;

        // get the most likely option
        int expMult = getExpMult();

        return P[expMult-firstMult];
}

double Multiplicity::operator[](int mult) const
{
        if (mult < firstMult)
                return 0.0;
        if (mult >= firstMult + (int)P.size())
                return 0.0;
        return exp(P[mult-firstMult]);
}

void Multiplicity::normalize()
{
        if (P.empty())
                return;

        double logSum = P.front();
        for (size_t i = 1; i < P.size(); i++) {
                double m = max<double>(logSum, P[i]);
                logSum = m + log( exp(P[i]-m) + exp(logSum-m) );
        }

        for (size_t i = 0; i < P.size(); i++)
                P[i] = P[i] - logSum;
}

std::ostream &operator<<(std::ostream &out, const Multiplicity &m)
{
        for (int i = 0; i < (int)m.P.size(); i++)
                out << i+m.firstMult << "\t" << m.P[i] << "\t";

        return out;
}

// ============================================================================
// COVERAGE MODEL
// ============================================================================

CovModel::CovModel(double errorLambda, double errorODF, double lambda,
                   double ODF, const std::vector<double>& w) :
        errorLambda(errorLambda), errorODF(errorODF), lambda(lambda), ODF(ODF)
{
        assert(w.size() > 2); // we want at least a weight for errors, mult 1 and mult > 1
        // store the logarithm of the weights
        logw.resize(w.size());
        const double pseudoCount = 1.0;
        for (size_t i = 0; i < logw.size(); i++)
                logw[i] = log(max<double>(w[i], pseudoCount));
}

CovModel::CovModel(const std::string& filename)
{
        read(filename);
}

int CovModel::getExpMult(double obsCov) const
{
        // Note that this routine will return either mult or mult+1. Depending
        // on the weights and lambda values, this might be an approximation.
        int mult = int(obsCov / lambda);

        double lambdaLo = (mult == 0) ? errorLambda : mult * lambda;
        double lambdaHi = (mult+1) * lambda;

        double wLo = (mult < (int)logw.size()) ? logw[mult] : logw.back();
        double wHi = (mult+1 < (int)logw.size()) ? logw[mult+1] : logw.back();

        double ODFLo = (mult == 0) ? errorODF : ODF;

        if ((wLo + Util::logNegbinomialPDF(obsCov, lambdaLo, ODFLo*lambdaLo)) <
            (wHi + Util::logNegbinomialPDF(obsCov, lambdaHi, ODF*lambdaHi)))
                mult++;

        return mult;
}

Multiplicity CovModel::getMultSoft(double obsCov, int numAlt) const
{
        int center = getExpMult(obsCov);
        int lo = max<int>(center - numAlt, 0);
        int hi = center + numAlt;

        vector<double> P(hi - lo + 1);
        for (int mult = lo; mult <= hi; mult++)
                P[mult-lo] = getLogProb(obsCov, mult);

        return Multiplicity(lo, P);
}

double CovModel::getLogProb(double obsCov, int mult) const
{
        double myLambda = (mult == 0) ? errorLambda : mult * lambda;
        double myODF = (mult == 0) ? errorODF : ODF;
        double w = (mult < (int)logw.size()) ? logw[mult] : logw.back();
        return w + Util::logNegbinomialPDF(obsCov, myLambda, myODF * myLambda);
        //return log(exp(w + Util::logNegbinomialPDF(obsCov, myLambda, myODF * myLambda)) + 1e-3);
}

void CovModel::write(const std::string& covFilename) const
{
        ofstream covFile(covFilename.c_str());

        covFile << errorLambda << "\t" << errorODF << "\t"
                << lambda << "\t" << ODF << "\t" << logw.size() << "\n";

        for(size_t i = 0; i < logw.size(); i++) {
                covFile << exp(logw[i]);
                covFile << (i+1 < logw.size() ? "\t" : "\n");
        }

        covFile.close();
}

void CovModel::read(const std::string& covFilename)
{
        ifstream ifs(covFilename.c_str());
        ifs >> errorLambda >> errorODF >> lambda >> ODF;
        size_t w_size;
        ifs >> w_size;
        logw.resize(w_size);
        for(size_t i = 0; i< w_size; i++) {
                double w;
                ifs >> w;
                logw[i] = log(w);
        }
}

void CovModel::writeGnuplot(const string& baseFilename,
                            const map<int, double>& hist) const
{
        ofstream ofs((baseFilename + ".dat").c_str());

        for (int c = 1; c < (logw.size() + 1) * lambda; c++) {
                ofs << c << "\t";
                for (int m = 0; m < logw.size(); m++)
                        ofs << exp(getLogProb(c, m)) << "\t";
                ofs << (hist.find(c) == hist.end() ? 0 : hist.at(c));
                ofs << "\n";
        }

        ofs.close();

        ofs.open((baseFilename + ".gnuplot").c_str());
        ofs << "set terminal pdf\n"
            << "set output \"" << baseFilename << "spectrum.pdf\"\n"
            << "set title \"Coverage histogram with fitted model\"\n"
            << "set xlabel \"coverage\"\n"
            << "set ylabel \"number of observations\"\n"
            << "set xrange [0:" << int(logw.size() * lambda) << "]\n"
            << "set yrange [0:" << int(exp(getLogProb(lambda, 1))*2.0) << "]\n"
            << "plot \"" << baseFilename << ".dat\" using 1:2 title \'seq. errors\' with lines,\\\n";
        string st = "($2";
        for (int m = 1; m < logw.size(); m++){
                st += ("+$" + to_string(m+2));
                ofs << "\t\"" << baseFilename << ".dat\" using 1:" << m+2 << " title \'mult = " << m << "\' with lines,\\\n";
        }
        st += ")";
        ofs << "\t\"" << baseFilename << ".dat\" using 1:" << logw.size()+2 << " title \'k-mers\' with boxes,\\\n";
        ofs << "\t\"" << baseFilename << ".dat\" using 1:" << st << " title \'mixture\' with lines\n";
        //ofs << "pause -1\n";
}

std::ostream &operator<<(std::ostream &out, const CovModel &c)
{
        out << fixed;
        out.precision(3);
        out << "[lambda: " << c.errorLambda << " (" << c.errorODF << ") - "
            << c.lambda << " (" << c.ODF << ")]\n\tweights:";
        out.precision(1);
        for (size_t i = 0; i < c.logw.size(); i++)
                out << " " << exp(c.logw[i]) << " (" << i << ")";
        out << "]";

        return out;
}

// ============================================================================
// DBGRAPH
// ============================================================================

double DBGraph::getInitialKmerCovEstimate(double errLambda, double p) const
{
        // sanity checks
        assert(errLambda > 0.0);
        assert(p > 0.0);
        assert(p < 1.0);

        // Given a Poisson distribution for the error model, find a cutoff
        // value for the coverage for which the probability of observing
        // a coverage is less than p under this error model
        double cutoff = ceil(errLambda);
        for ( ; cutoff < 10.0 * errLambda; cutoff++)
                if (Util::poissonPDF((unsigned int)cutoff, errLambda) < p)
                        break;

        double totCoverage = 0, totSize = 0;
        for ( NodeID id = 1; id <= numNodes; id++ ) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                if (node.getAvgCov() < cutoff)
                        continue;
                totCoverage += node.getCov();
                totSize += node.getMarginalLength();
        }

        if (totSize > 0)
                return (double)totCoverage / (double)totSize;
        else
                return 2.0 * errLambda;      // pathological case
}

void DBGraph::covCount(const FastQRecord& rr, const KmerNPPTable& table)
{
        const string& read = rr.getRead();
        const string& qual = rr.getQual();

        NodePosPair prevNpp; NodePosition prevOff = 0;
        for (KmerIt it(read); it.isValid(); it++) {
                NodePosPair npp = table.find(it.getKmer());
                if (!npp.isValid())
                        continue;

                // increase the node coverage
                NodeID id = npp.getNodeID();
                SSNode n = getSSNode(id);
                double score = Util::phred2prob(qual, it.getOffset(),
                                                it.getOffset() + Kmer::getK());
                n.incCov(score);
#ifndef AVG_COV
		n.incCount(npp.getPosition(), 1);
#endif

                // if current kmer succeeds a valid previous kmer
                if ((prevOff+1 == it.getOffset()) && (crossesArc(prevNpp, npp))) {
                        score = Util::phred2prob(qual, prevOff,
                                                 prevOff + Kmer::getK() + 1);
                        getSSNode(prevNpp.getNodeID()).rightArc(id)->incCov(score);
                        // palindromic arcs exist only once! Don't add coverage!
                        if (prevNpp.getNodeID() != - id)
                                n.leftArc(prevNpp.getNodeID())->incCov(score);
#ifndef AVG_COV
		} else if (prevOff+1 == it.getOffset() && prevNpp.getNodeID() == id && it.getOffset() != 0) {
		  n.incArcCount(npp.getPosition()-1, 1);
#endif
                }

                prevNpp = npp;
                prevOff = it.getOffset();
        }
}

void DBGraph::covCountThread(FastQReader& inputs,
                             const KmerNPPTable& table)
{
        // local storage of reads
        vector<FastQRecord> readBuf;

        size_t chunkID;
        while (inputs.getNextChunk(readBuf, chunkID))
                for (const auto& readRecord : readBuf)
                        covCount(readRecord, table);
}

void DBGraph::getCovFromReads(LibraryContainer &inputs, const KmerNPPTable& table)
{
        // reset all coverages (might have been loaded from BCALM)
        for (size_t i = 1; i <= numNodes; i++) {
                getSSNode(i).setCov(0);
#ifndef AVG_COV
		getSSNode(i).clearCounts();
		getSSNode(i).clearArcCounts();
#endif
	}
	
        for (size_t i = 1; i <= numArcs; i++) {
                arcs[i].setCov(0);
	}

        // initialize Phred score lookup table
        if (settings.useQual()) {
                cout << "Using Phred quality scores (ASCII base="
                     << settings.getPhredBase() << ") to weigh k-mers\n";
                Util::enablePhred(settings.getPhredBase());
        } else {
                cout << "Not using Phred quality scores\n";
        }

        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        const size_t ws = settings.getThreadIOWorkSize();
        for (size_t i = 0; i < inputs.size(); i++) {
                string fn1, fn2;
                tie(fn1, fn2) = inputs.getFilename(i);

                FastQReader myReader(fn1, fn2);
                myReader.startReaderThread(ws, ws * settings.getNumThreads());

                // start worker threads
                vector<thread> workerThreads(numThreads);
                for (size_t i = 0; i < workerThreads.size(); i++)
                        workerThreads[i] = thread(&DBGraph::covCountThread,
                                                  this, ref(myReader), cref(table));

                // wait for worker threads to finish
                for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

                myReader.joinReaderThread();
        }

        size_t numZero = 0;
        for (ArcID id = 1; id <= numArcs; id++)
                if (arcs[id].getCov() == 0)
                        numZero++;

        if (numZero > 0)
                cerr << "WARNING: found " << numZero << " arcs with coverage 0\n";
}
