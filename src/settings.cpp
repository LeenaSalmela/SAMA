/******************************************************************************
 *   Copyright (C) 2024 Leena Salmela (leena.salmela@helsinki.fi)             *
 *   This file has been modified for SAMA                                     *
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

#include "settings.h"
#include "util.h"

#include <thread>
#include <iostream>
#include <fstream>

using namespace std;

void Settings::printProgramVersion() const
{
        cout << "SAMA, a de Bruijn graph assembler with guaranteed correctness -- version "
             << SAMA_MAJOR_VERSION << "." << SAMA_MINOR_VERSION
             << "." << SAMA_PATCH_LEVEL << "\n";

        cout << "Copyright (C) 2019-2022 Jan Fostier and Aranka Steyaert\n"
                "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n"
                "This is free software; see the source for copying conditions. "
                "There is NO\nwarranty; not even for MERCHANTABILITY or "
                "FITNESS FOR A PARTICULAR PURPOSE." << endl;
}

void Settings::printUsage() const
{
        cout <<

        "Usage: sama [options] unitigs.fa readfile\n\n  File "
        "\"unitigs.fa\" is the input de Bruijn graph produced by BCALM 2\n"
        "  File \"readfile\" contains a list of input FASTQ files\n\n"

        " [options without arg]\n"
        "  -use-qual\t\tUse Phred quality scores to weigh k-mer coverage. When using this\n"
        "            \t\toption, it is advised to also pass the -abundance-min parameter\n"
        "  -approx-inf\t\tUse loopy belief propagation on the whole de Bruijn graph whenever\n"
        "              \t\tthe CRF is used for computation\n"
        "  -map-assignment\t\t Use in combination with approx-inf; compute a MAP assignment\n\n"

        "  -help\t\t\tdisplay help page\n"
        "  -version\t\tdisplay version\n\n"

        " [options with 1 arg]\n"
        "  -num-threads\t\tNumber of threads [default = #cores]\n\n"

        "  -abundance-min\tMin abundance threshold as passed to BCALM 2 [default = auto]\n\n"

        "  -crf-nb-size\t\tCRF neighborhood size [default = 5]\n"
        "  -crf-margin\t\tCRF factor margin [default = 2]\n"
        "  -crf-flow\t\tCRF flow conservation strength [default = 1.0e7]\n"
        "  -crf-max-fact\t\tCRF maximum factor size [default = 1.0e6]\n\n"

        "  -mm-coverage\t\tInitial coverage est. in mixture model [default = auto]\n"
        "  -mm-err-cov\t\tInitial error coverage est. in mixture model [default = 1.0]\n"
        "  -mm-odf\t\tInitial overdispersion factor est. in mixture model [default = 1.5]\n"
        "  -mm-components\tNumber of components in mixture model [default = 6]\n\n"

        "  -em-max-iter\t\tMaximum number of EM iterations [default = 25]\n"
        "  -em-conv-eps\t\tRelative EM convergence epsilon [default = 1.0e-3]\n"
        "  -em-train-size\tNumber of nodes/arcs to use for EM training [default = 1e4]\n\n"

        "  -phred-base\t\tASCII value corresponding to Phred score Q = 0 [default = 33]\n\n"
        "  -vis-subgraph\t\tCentral nodeID around which a subgraph will be visualised [default = 0]\n"
        "               \t\t if value = 0, no visualisation and normal pipeline will be run\n"

	"  -misassembly-likelihood\tLikelihood of a misassembly [default=1.0e-6]\n\n"

        "Report bugs to Leena Salmela <leena.salmela@helsinki.fi>" << endl;
}

Settings::Settings(int argc, char ** argv) : useQualFlag(false), approxInfFlag(false),
        mapAssignment(false), singleCRF(false), numThreads(thread::hardware_concurrency()),
        abundanceMin(-1), crfDepth(5), crfMargin(2), crfFlow(1e7), crfMaxFact(1e6),
        mmCov(-1.0), mmErrCov(1.0), mmODF(1.5), mmComponents(6), emMaxIter(25),
	emConvEps(1e-3), emTrainSize(1e4), phredBase(33), visGraphNode(0), misassLh(1e-10)
{
        const int reqArguments = 2;     // not counting argument 0

        // no arguments provided
        if (argc <= 1) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // only one argument provided
        if (argc == 2) {
                string arg(argv[1]);

                if (arg == "-help") {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if (arg == "-version") {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        // process optional arguments (we know there are at least two arguments)
        for (int i = 1; i < argc-reqArguments; i++) {
                string arg(argv[i]);

                // process options without arguments
                if (arg == "-use-qual") {
                        useQualFlag = true;
                        continue;
                } else if (arg == "-approx-inf") {
                        approxInfFlag = true;
                        continue;
                } else if (arg == "-map-assignment") {
                        mapAssignment = true;
                        continue;
                } else if (arg == "-single-crf"){
                        singleCRF = true;
                        continue;
                } else if (arg == "-help") {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if (arg == "-version") {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                } else if (i == argc-reqArguments-1) {
                        cerr << "Unknown option or missing argument: "
                             << argv[i] << endl;
                        printUsage();
                        exit(EXIT_FAILURE);
                }

                // process options with arguments
                if (arg == "-num-threads") {
                        numThreads = atoi(argv[i+1]);
                } else if (arg == "-abundance-min") {
                        abundanceMin = atoi(argv[i+1]);
                } else if (arg == "-crf-nb-size") {
                        crfDepth = atoi(argv[i+1]);
                } else if (arg == "-crf-margin") {
                        crfMargin = atoi(argv[i+1]);
                } else if (arg == "-crf-flow") {
                        crfFlow = atof(argv[i+1]);
                } else if (arg == "-crf-max-fact") {
                        crfMaxFact = atof(argv[i+1]);
                } else if (arg == "-mm-coverage") {
                        mmCov = atof(argv[i+1]);
                } else if (arg == "-mm-err-cov") {
                        mmErrCov = atof(argv[i+1]);
                } else if (arg == "-mm-odf") {
                        mmODF = atof(argv[i+1]);
                } else if (arg == "-mm-components") {
                        mmComponents = atoi(argv[i+1]);
                } else if (arg == "-em-max-iter") {
                        emMaxIter = atoi(argv[i+1]);
                } else if (arg == "-em-conv-eps") {
                        emConvEps = atof(argv[i+1]);
                } else if (arg == "-em-train-size") {
                        emTrainSize = atof(argv[i+1]);
                } else if (arg == "-phred-base") {
                        phredBase = atoi(argv[i+1]);
                } else if (arg == "-vis-subgraph") {
                        visGraphNode = atoi(argv[i+1]);
                } else if (arg == "-misassembly-likelihood") {
                        misassLh = atof(argv[i+1]);
                } else {
                        cerr << "Unknown argument: " << argv[i] << endl;
                        printUsage();
                        exit(EXIT_FAILURE);
                }

                i++;    // if we reach this point, an argument was processed
        }
        
        graphFilename = argv[argc-2];
        readFilename = argv[argc-1];

        // check graph file
        if (graphFilename.empty()) {
                cerr << "Specify input de Bruijn graph file" << endl;
                printUsage();
                exit(EXIT_FAILURE);
        }

        if (!Util::fileExists(graphFilename))
                throw runtime_error("cannot open file " + graphFilename);

        // check read file
        if (readFilename.empty()) {
                cerr << "Specify input read file" << endl;
                printUsage();
                exit(EXIT_FAILURE);
        }

        if (!Util::fileExists(readFilename))
                throw runtime_error("cannot open file " + readFilename);
        
        // perform sanity checks on input paramters
        if (numThreads < 1)
                throw runtime_error("Number of threads must be >= 1");
        if (crfDepth < 0 || crfDepth > 10)
                if (! approxInfFlag)
                        throw runtime_error("CRF subgraph depth must be "
                                            "in range [0..10]");
        if (mapAssignment && ! approxInfFlag)
                throw runtime_error("MAP assignments not supported with exact inference on neightbourhoods");
        if (crfMargin < 0 || crfMargin > 5)
                throw runtime_error("CRF factor margin must be "
                                    "in range [0..5]");
        if (crfFlow < 1.0)
                throw runtime_error("CRF flow conservation strength "
                                    "must be >= 1.0");
        if (crfMaxFact < 1e3)
                throw runtime_error("CRF maximum factor size must "
                                    "be >= 1e3");
        // Note: negative avgCov means "auto detect" and is hence allowed
        if (mmCov >= 0.0 && mmCov < 1.0)
                throw runtime_error("Mixture model initial coverage "
                                    "must be >= 1.0");
        if (mmCov >= 0.0 && mmCov <= mmErrCov)
                throw runtime_error("Mixture model initial coverage must be "
                                    "higher than the initial error coverage");
        if (mmErrCov <= 0.0)
                throw runtime_error("Mean sequencing error coverage "
                                    "must be > 0.0");
        if (mmODF < 1.0)
                throw runtime_error("Mixture model initial overdispersion "
                                    "factor must be >= 1.0");
        if (mmComponents < 2 || mmComponents > 15)
                throw runtime_error("Mixture model number of components must "
                                    "in range [2..15]");
        if (emMaxIter < 1)
                throw runtime_error("Maximum number of EM iterations "
                                    "must be >= 1");
        if (emConvEps >= 1.0 || emConvEps < 1e-8)
                throw runtime_error("EM relative converage tolerance "
                                    "should be in range [1e-8...1[");
        if (emTrainSize < 1e2)
                throw runtime_error("EM traing size should be >= 1e2");
        if (phredBase < 0 || phredBase > 255)
                throw runtime_error("ASCII value for Phred score Q = 0 must "
                                    "be in range [0..255]");
        if (useQualFlag && abundanceMin == -1) {
                cerr << "WARNING: no -abundance-min argument passed while using q-mers.\n"
                     << "SAMA cannot automatically infer the correct value.\n"
                     << "Please consider providing the appropriate -abudance-min argument to SAMA.\n";
        }

	if (misassLh <= 0 || misassLh >= 1.0) {
	  throw runtime_error("Misassembly likelihood must be larger than 0.0 and smaller than 1.0");
	}
}
