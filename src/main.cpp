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

#include <cstdlib>
#include <numeric>

#include "dbgraph.h"
#include "settings.h"
#include "refcomp.h"
#include "crfmult.h"
#include "readaln.h"

using namespace std;

void populateNodeMult(NodeMap<Multiplicity>& nodeMult,
                      const vector<NodeRep>& nodes)
{
        nodeMult = NodeMap<Multiplicity>(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++)
                nodeMult[nodes[i]] = Multiplicity();
}

void populateEdgeMult(EdgeMap<Multiplicity>& edgeMult,
                      const vector<EdgeRep>& edges)
{
        edgeMult = EdgeMap<Multiplicity>(edges.size());
        for (size_t i = 0; i < edges.size(); i++)
                edgeMult[edges[i]] = Multiplicity();
}

void trainModel(DBGraph& dBG, Settings& settings,
                CovModel& nodeModel, CovModel& edgeModel, bool fixedZero = false)
{
        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());

        NodeMap<Multiplicity> nodeMult;
        EdgeMap<Multiplicity> edgeMult;
        size_t trainSize = settings.getEMTrainingSize();
        populateNodeMult(nodeMult, dBG.getNodeReps(trainSize));
        populateEdgeMult(edgeMult, dBG.getEdgeReps(trainSize));

        int nIters = myCRFMult.computeMultEM(nodeMult, nodeModel,
                                             edgeMult, edgeModel,
                                             settings.getEMConvEps(),
                                             settings.getEMMaxIter(),
                                             false, false, fixedZero); //no approximate inference during model training

        if (nIters <= settings.getEMMaxIter())
                cout << "EM algorithm converged after "
                     << nIters << " iterations\n";
        else
                cout << "WARNING: maximum number of iterations reached. "
                     << "Convergence is not guaranteed." << endl;

        map<int, double> nodeHist;
        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                double f = node.getAvgCov() - floor(node.getAvgCov());
                nodeHist[node.getAvgCov()] += (1.0-f);
                nodeHist[node.getAvgCov() + 1] += f;
        }

        string app = fixedZero? ".st3" : ".st2";
        nodeModel.writeGnuplot("histogram.node" + app, nodeHist);

        map<int, double> edgeHist;
        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(dBG.getArcID(it.first));
                double f = arc.getCov() - floor(arc.getCov());
                edgeHist[arc.getCov()] += (1.0-f);
                edgeHist[arc.getCov() + 1] += f;
        }

        edgeModel.writeGnuplot("histogram.edge" + app, edgeHist);
}

void computeMultiplicities(DBGraph& dBG, Settings& settings,
                           NodeMap<Multiplicity>& nodeMult,
                           EdgeMap<Multiplicity>& edgeMult,
                           bool approx_inf=false, bool map=false)
{
        // get the node multiplicities
        CovModel nodeModel(settings.getNodeModelFilename());
        //nodeModel.setWeightUniform(1.0);
        CovModel edgeModel(settings.getEdgeModelFilename());
        //edgeModel.setWeightUniform(1.0);

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());


        populateNodeMult(nodeMult, dBG.getNodeReps(dBG.getNumValidNodes()));
        populateEdgeMult(edgeMult, dBG.getEdgeReps(dBG.getNumValidArcs()));

        if(approx_inf)
                myCRFMult.approxMultAll(nodeMult, nodeModel, edgeMult, edgeModel,
                                        map, settings.useSingleCRF());
        else
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);
}

void computeMultiplicities(DBGraph& dBG, Settings& settings,
                           NodeMap<Multiplicity>& nodeMult,
                           EdgeMap<Multiplicity>& edgeMult,
                           const CovModel& nodeModel, const CovModel& edgeModel,
                           bool approx_inf=false, bool map=false)
{
        
        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());
        
        
        populateNodeMult(nodeMult, dBG.getNodeReps(dBG.getNumValidNodes()));
        populateEdgeMult(edgeMult, dBG.getEdgeReps(dBG.getNumValidArcs()));
        
        if(approx_inf)
                myCRFMult.approxMultAll(nodeMult, nodeModel, edgeMult, edgeModel,
                                        map, settings.useSingleCRF());
        else
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);
}

void readTrueMultiplicities(NodeMap<int>& trueNodeMult,
                            EdgeMap<int>& trueEdgeMult,
                            string tnmF, string temF){
        ifstream ifs(tnmF);
        if (!ifs)
                cerr << "Could not find true nodes multiplicities file...\n"
                "True node multiplicities will be set to -1\n";
        
        while (ifs) {
                NodeID nodeID; int m;
                ifs >> nodeID >> m;
                if (!ifs)
                        break;
                auto it = trueNodeMult.find(nodeID);
                if (it != trueNodeMult.end())
                        it->second = m;
        }
        ifs.close();
        
        ifs.open(temF);
        if (!ifs)
                cerr << "Could not find true edges multiplicities file...\n"
                "True edge multiplicities will be set to -1\n";
        
        while (ifs) {
                NodeID srcID, dstID; int m;
                ifs >> srcID >> dstID >> m;
                if (!ifs)
                        break;
                auto it = trueEdgeMult.find(EdgeRep(srcID, dstID));
                if (it != trueEdgeMult.end())
                        it->second = m;
        }
        ifs.close();
}


void drawFullCytoGraph(Settings& settings, DBGraph& dBG) 
{
        vector<NodeID> nID; vector<EdgeID> eID;
        dBG.getGraph(nID, eID);
        
        vector<NodeRep> nodes;
        vector<EdgeRep> edges;
        NodeMap<Multiplicity> estNodeMult;
        EdgeMap<Multiplicity> estEdgeMult;
        
        ifstream ifs(settings.getEstNodeMultFilename());
        if (!ifs){
                cerr << "Could not find estimated nodes multiplicities file...\n"
                "Cannot write cytoscape graph\n";
                return;
        }
        
        while (ifs) {
                NodeID nodeID; int m; double lp; double d1; size_t d2;
                ifs >> nodeID >> m >> lp >> d1 >> d2;
                if (!ifs)
                        break;
                nodes.push_back(NodeRep(nodeID));
                estNodeMult[NodeRep(nodeID)] = Multiplicity(m, {lp,log(1.0-exp(lp))});
        }
        ifs.close();
        
        ifs.open(settings.getEstEdgeMultFilename());
        if (!ifs){
                cerr << "Could not find estimated edge multiplicities file...\n"
                "Cannot write cytoscape graph\n";
                return;
        }
        
        while (ifs) {
                NodeID srcID; NodeID dstID; int m; double lp; double d1;
                ifs >> srcID >> dstID >> m >> lp >> d1;
                if (!ifs)
                        break;
                edges.push_back(EdgeRep(srcID, dstID));
                estEdgeMult[EdgeRep(srcID, dstID)] = Multiplicity(m, {lp, log(1.0-exp(lp))});
        }
        ifs.close();
        
        // If they exist, also load the true multiplicities from disk
        // Note: for efficiency, story only true multiplicities that are needed
        NodeMap<int> trueNodeMult;
        transform(nodes.begin(), nodes.end(),
                  inserter(trueNodeMult, trueNodeMult.end()),
                  [](const NodeRep& nr) { return make_pair(nr, Multiplicity(-1)); });
        
        EdgeMap<int> trueEdgeMult;
        transform(edges.begin(), edges.end(),
                  inserter(trueEdgeMult, trueEdgeMult.end()),
                  [](const EdgeRep& er) { return make_pair(er, Multiplicity(-1)); });
        
        ifs.open("truemult.node");
        if (!ifs)
                cerr << "Could not find true nodes multiplicities file...\n"
                "True node multiplicities will be set to -1 in "
                "resulting Cytoscape graph\n";
        
        while (ifs) {
                NodeID nodeID; int m;
                ifs >> nodeID >> m;
                if (!ifs)
                        break;
                auto it = trueNodeMult.find(nodeID);
                if (it != trueNodeMult.end())
                        it->second = m;
        }
        ifs.close();
        
        ifs.open("truemult.edge");
        if (!ifs)
                cerr << "Could not find true edges multiplicities file...\n"
                "True edge multiplicities will be set to -1 in "
                "resulting Cytoscape graph\n";
        
        while (ifs) {
                NodeID srcID, dstID; int m;
                ifs >> srcID >> dstID >> m;
                if (!ifs)
                        break;
                auto it = trueEdgeMult.find(EdgeRep(srcID, dstID));
                if (it != trueEdgeMult.end())
                        it->second = m;
        }
        ifs.close();
        
        cout << "Writing Cytoscape graph... " << endl;
        dBG.writeCytoscapeGraph("Cytograph.full",
                                nID, eID,
                                estNodeMult, estEdgeMult,
                                trueNodeMult, trueEdgeMult);
}

void visualiseSubgraph(Settings& settings){
        cout << "\nVisualising neighbourhood of size " 
                << settings.getCRFDepth() << " around node "
                << settings.getVisGraphNode() 
                << "\n==================================================\n" << endl;
        
        NodeID noi = settings.getVisGraphNode();
        string ifn;
        bool corrected;
        if ( ! settings.stageThreeNecessary() ){
                cout << "\nUsing cleaned de Bruijn graph of stage 3\n" << endl;
                ifn = settings.getStage3GraphFilename();
                corrected  = true;
        } else {
                cout << "\nUsing uncleaned de Bruijn graph of stage 1\n" << endl;
                ifn = settings.getStage1GraphFilename();
                corrected = false;
        }
        
        DBGraph dBG(settings);
        dBG.loadBinary(ifn);
        
        if((noi < -dBG.getNumNodes()) || (noi > dBG.getNumNodes()))
                throw runtime_error("Specified node id of " + to_string(noi) +
                                        " does not exist in the graph");
        
        if(! dBG.getSSNode(noi).isValid())
                throw runtime_error("Specified node id of " + to_string(noi) +
                                        " does not exist in the graph");
        
        CovModel nodeModel(settings.getNodeModelFilename(corrected));
        CovModel edgeModel(settings.getEdgeModelFilename(corrected));
        
        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;
        
        vector<NodeID> nodes;
        vector<pair<NodeID, NodeID> > edges;
        
        dBG.getSubgraph(noi, nodes, edges, settings.getCRFDepth());
        
        vector<NodeRep> nodeReps(nodes.begin(), nodes.end());
        vector<EdgeRep> edgeReps(edges.begin(), edges.end());
        
        // Convert estimated node/edge multiplicities to lookup table
        NodeMap<Multiplicity> estNodeMult;
        EdgeMap<Multiplicity> estEdgeMult;
        
        if (settings.approxInf()){
                CRFSolver myCRFSolver(dBG, settings.getCRFMargin(),
                                      settings.getCRFMaxFactSize(),
                                      settings.getCRFFlowStrength());
                myCRFSolver.approxSubgraphMult(noi,
                                               estNodeMult, estEdgeMult,
                                               nodeModel, edgeModel, settings.getCRFDepth(), settings.computeMAP());
        }else{
                CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                                  settings.getCRFMargin(),
                                  settings.getCRFMaxFactSize(),
                                  settings.getCRFFlowStrength(),
                                  settings.getNumThreads(),
                                  settings.getThreadGraphWorkSize());
                populateNodeMult(estNodeMult, nodeReps);
                populateEdgeMult(estEdgeMult, edgeReps);
                
                myCRFMult.computeMult(estNodeMult, nodeModel, estEdgeMult, edgeModel);
        }
        
        NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;
        transform(nodeReps.begin(), nodeReps.end(),
                  inserter(trueNodeMult, trueNodeMult.end()),
                  [](const NodeRep& nr) { return make_pair(nr, Multiplicity(-1)); });
        transform(edgeReps.begin(), edgeReps.end(),
                  inserter(trueEdgeMult, trueEdgeMult.end()),
                  [](const EdgeRep& er) { return make_pair(er, Multiplicity(-1)); });
        
        readTrueMultiplicities(trueNodeMult, trueEdgeMult, 
                               settings.getTrueNodeMultFilename() + (corrected ? "" : ".st1"),
                               settings.getTrueEdgeMultFilename() + (corrected ? "" : ".st1"));
        
        cout << "Writing Cytoscape graph... " << endl;
        dBG.writeCytoscapeGraph("cytgraph" + to_string(noi) + "nb" + to_string(settings.getCRFDepth()),
                                nodes, edges,
                                estNodeMult, estEdgeMult,
                                trueNodeMult, trueEdgeMult);
}


void stageOne(Settings& settings, LibraryContainer& libraries)
{
        cout << "\nEntering stage 1\n";
        cout << "================\n" << endl;

        if (!settings.stageOneNecessary()) {
                cout << "File " << settings.getStage1GraphFilename()
                     << " exists. Skipping stage 1...\n";
                return;
        }

        Util::startChrono();

        DBGraph dBG(settings);

        // load the graph from BCALM2 file
        dBG.loadBCalm(settings.getGraphFilename());

        // create a <kmer, nodePosPair> table
        KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);

        // stream the reads to get the coverage
        dBG.getCovFromReads(libraries, table);
        
        // write stage 1 binary output file
        dBG.writeBinary(settings.getStage1GraphFilename());

        if (Util::fileExists("genome.fasta")) {
                //dBG.sanityCheck();

                RefComp refComp(dBG, settings, table, "genome.fasta");
                //refComp.writeAlignedSeqs("theoretical_contigs.fasta");
		std::vector<std::vector<AlnSegment> > aln {};
                refComp.alignSequences(aln);

                NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;
                refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);

                vector<NodeRep> nodes = dBG.getNodeReps(dBG.getNumValidNodes());
                vector<EdgeRep> edges = dBG.getEdgeReps(dBG.getNumValidArcs());

                ofstream nodeOFS(settings.getTrueNodeMultFilename()+".st1");
                ofstream edgeOFS(settings.getTrueEdgeMultFilename()+".st1");

                // write the true node multiplicities
                 for (const auto& it : trueNodeMult)
                        nodeOFS << it.first << "\t" << it.second << "\n";
                 nodeOFS.close();
                
                // write the true edge multiplicities
                 for (const auto& it : trueEdgeMult)
                        edgeOFS << it.first.getSrcID() << "\t"
                                << it.first.getDstID() << "\t"
                                << it.second << "\n";
                 edgeOFS.close();
        }

        cout << "Stage 1 finished in " << Util::stopChronoStr() << endl;
}

void stageTwo(Settings& settings)
{
        cout << "\nEntering stage 2\n";
        cout << "================\n" << endl;

        if (!settings.stageTwoNecessary()) {
                cout << "Files " << settings.getNodeModelFilename()
                     << " and " << settings.getEdgeModelFilename()
                     << " exist. Skipping stage 2...\n";
		return;
        }

        Util::startChrono();

        // load stage 1 de Bruijn graph
        DBGraph dBG(settings);
        dBG.loadBinary(settings.getStage1GraphFilename());

        // TODO: uniform initialization of weights is best ??
        vector<double> wN(settings.getNumComp(), 1.0);
        vector<double> wE(settings.getNumComp(), 1.0);

        double initErrCov = settings.getInitErrCov();
        double initCov = settings.getInitCov();
        double initODF = settings.getInitODF();

        if (initCov < 0.0)
                initCov = dBG.getInitialKmerCovEstimate(initErrCov, 0.01);
        cout << "Initial coverage estimate: " << fixed << initCov << endl;

        CovModel edgeModel(initErrCov, initODF, initCov, initODF, wE);
        CovModel nodeModel(initErrCov, initODF, initCov, initODF, wN);

        // train and write coverage models
        // no approximate inference here?
        trainModel(dBG, settings, nodeModel, edgeModel);
        
        nodeModel.write(settings.getNodeModelFilename());
        edgeModel.write(settings.getEdgeModelFilename());

        cout << "Stage 2 finished in " << Util::stopChronoStr() << endl;
}

void stageAssemble(Settings& settings)
{
  if (!Util::fileExists(settings.getThresholdFilename())) {
    cout << "File " << settings.getThresholdFilename()
	 << " does not exist. Compute the thresholds file and rerun.\n";
    return;
  }

  // Read the thresholds

  // Examine each arc and remove it if it does not meet the threshold

  // Generate unitigs
}

void stageThreeCorrectGraph(Settings& settings)
{
        cout << "\nEntering stage 3\n";
        cout << "================\n" << endl;

        if (Util::fileExists(settings.getStage3GraphFilename())) {
                cout << "File " << settings.getStage3GraphFilename()
                     << " exist. Skipping stage 3...\n";
		return;
        }

        Util::startChrono();

        // load stage 1 graph
        DBGraph dBG(settings);
        string ifn = settings.getStage1GraphFilename();
        cout << "Loading graph from file: " << ifn << "..."; cout.flush();
        Util::startChrono();
        dBG.loadBinary(ifn);
        cout << "\n\tLoaded " << dBG.getNumNodes() << " nodes and "
             << dBG.getNumArcs() << " arcs (" << Util::stopChronoStr() << ")\n";

        // load node/edge coverage model
        CovModel nodeModel(settings.getNodeModelFilename());
        CovModel edgeModel(settings.getEdgeModelFilename());

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getCRFWorkSize());

        // remove arcs with zero coverage
        cout << "Removing nodes/arcs with zero coverage..." << endl;
        dBG.removeCoverage(0.0, 0);
        dBG.concatenateNodes();

        cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
             << dBG.getNumValidArcs() << " arcs" << endl;

        // remove low-coverage nodes if the FLOW is OK
        const int numRounds = 4;
        for (int i = 1; i <= numRounds; i++)
        {
                double fraction = double(i) / double(numRounds);
                double nodeCutoff = fraction * nodeModel.getCovCutOff(1);

                // a) TIPS
                {
                vector<NodeRep> nodeReps = dBG.getLowCovTips(nodeCutoff,
                                                             2*Kmer::getK());
                cout << "Selected " << nodeReps.size()
                     << " tips with coverage <= " << nodeCutoff << endl;

                vector<bool> flowOK(nodeReps.size());
                myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);

                size_t numRemove = count(flowOK.begin(), flowOK.end(), true);

                vector<NodeRep> toRemove;
                toRemove.reserve(numRemove);

                for (size_t i = 0; i < nodeReps.size(); i++)
                        if (flowOK[i])
                                toRemove.push_back(nodeReps[i]);

                cout << "\tRemoving " << numRemove << " nodes" << endl;
                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes()
                     << " nodes and " << dBG.getNumValidArcs() << " arcs\n";
                }

                // b) BUBBLES
                {
                vector<NodeRep> nodeReps = dBG.getLowCovBubbles(nodeCutoff);
                cout << "Selected " << nodeReps.size()
                     << " bubbles with coverage <= " << nodeCutoff << endl;

                vector<bool> flowOK(nodeReps.size());
                myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);

                size_t numRemove = count(flowOK.begin(), flowOK.end(), true);
                vector<NodeRep> toRemove;
                toRemove.clear();
                toRemove.reserve(numRemove);

                cout << "\tRemoving " << numRemove << " nodes" << endl;
                for (size_t i = 0; i < nodeReps.size(); i++)
                        if (flowOK[i])
                                toRemove.push_back(nodeReps[i]);

                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes()
                     << " nodes and " <<  dBG.getNumValidArcs() << " arcs\n";
                }

                // c) ALL
                {
                vector<NodeRep> nodeReps = dBG.getLowCovNodes(nodeCutoff);
                cout << "Selected " << nodeReps.size()
                     << " nodes with coverage <= " << nodeCutoff << endl;

                vector<bool> flowOK(nodeReps.size());
                myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);

                size_t numRemove = count(flowOK.begin(), flowOK.end(), true);
                vector<NodeRep> toRemove;
                toRemove.clear();
                toRemove.reserve(numRemove);

                cout << "\tRemoving " << numRemove << " nodes" << endl;
                for (size_t i = 0; i < nodeReps.size(); i++)
                        if (flowOK[i])
                                toRemove.push_back(nodeReps[i]);

                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes()
                     << " nodes and " <<  dBG.getNumValidArcs() << " arcs\n";
                }
        }
	{
        NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;
        if (Util::fileExists("genome.fasta")) {
                // create a <kmer, nodePosPair> table
                KmerNPPTable table(dBG, true);
                dBG.createKmerNPPTable(table);
                RefComp refComp(dBG, settings, table, "genome.fasta");
                //refComp.writeAlignedSeqs("theoretical_initial.fasta");
                refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);
        }
	}


        // remove low-coverage nodes and edges if the CRF is OK
        for (int i = 1; i <= numRounds+2; i++)
        {
                double fraction = double(i) / double(numRounds);
                double covCutoff = fraction * nodeModel.getCovCutOff(1);

                // a) nodes
                vector<NodeRep> nodeReps = dBG.getLowCovNodes(covCutoff);
                cout << "Selected " << nodeReps.size()
                     << " nodes with coverage <= " << covCutoff << endl;
                NodeMap<Multiplicity> nodeMult(nodeReps.size());
                for (size_t i = 0; i < nodeReps.size(); i++)
                        nodeMult[nodeReps[i]] = Multiplicity();
                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);

                vector<NodeRep> toRemove;
                for (const auto& it : nodeMult)
                        if (it.second.getExpMult() == 0)
                                toRemove.push_back(it.first);

                cout << "\tRemoving " << toRemove.size() << " nodes" << endl;
                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                // b) edges
                vector<EdgeRep> edgeReps = dBG.getLowCovEdges(covCutoff);
                cout << "Selected " << edgeReps.size()
                     << " edges with coverage <= " << covCutoff << endl;

                nodeMult.clear();
                edgeMult = EdgeMap<Multiplicity>(edgeReps.size());
                for (size_t i = 0; i < edgeReps.size(); i++)
                        edgeMult[edgeReps[i]] = Multiplicity();

                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);

                vector<EdgeRep> edgesToRemove;
                for (const auto& it : edgeMult)
                        if (it.second.getExpMult() == 0)
                                edgesToRemove.push_back(it.first);

                cout << "\tRemoving " << edgesToRemove.size() << " arcs" << endl;
                dBG.removeEdges(edgesToRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes()
                     << " nodes and " <<  dBG.getNumValidArcs() << " arcs\n";
        }

        // remove (high-coverage) tips if CRF is OK
        {
                vector<NodeRep> nodeReps = dBG.getLowCovTips(1e100,
                                                             2*Kmer::getK());
                cout << "Selected " << nodeReps.size() << " tips " << endl;

                NodeMap<Multiplicity> nodeMult(nodeReps.size());
                for (size_t i = 0; i < nodeReps.size(); i++)
                        nodeMult[nodeReps[i]] = Multiplicity();
                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel,
                                      edgeMult, edgeModel);

                vector<NodeRep> toRemove;
                for (size_t i = 0; i < nodeReps.size(); i++)
                        if (nodeMult[i].getExpMult() == 0)
                                toRemove.push_back(nodeReps[i]);

                cout << "\tRemoving " << toRemove.size() << " nodes" << endl;
                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes()
                     << " nodes and " <<  dBG.getNumValidArcs() << " arcs\n";
        }

	{
        NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;
        if (Util::fileExists("genome.fasta")) {
                // create a <kmer, nodePosPair> table
                KmerNPPTable table(dBG, true);
                dBG.createKmerNPPTable(table);
                RefComp refComp(dBG, settings, table, "genome.fasta");
                //refComp.writeAlignedSeqs("theoretical_initial.fasta");
                refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);
        }
	}

        // re-train and write coverage models (fix zero distribution)
        nodeModel.setZeroWeight(1.0);
        edgeModel.setZeroWeight(1.0);
        trainModel(dBG, settings, nodeModel, edgeModel, true);
        nodeModel.write(settings.getNodeModelFilename(true));
        edgeModel.write(settings.getEdgeModelFilename(true));

        dBG.defrag();
        dBG.sanityCheck();

        // write stage 3 binary output file
        dBG.writeBinary(settings.getStage3GraphFilename());
        cout << "Stage 3 finished in " << Util::stopChronoStr() << endl;
}

void StageFourComputeMult(Settings& settings)
{
        cout << "\nDetermining multiplicities all nodes and arcs in the de Bruijn graph\n";
        cout << "======================================================================\n" << endl;
        
        Util::startChrono();
        
        // Load the graph
        DBGraph dBG(settings);
        if ( ! Util::fileExists(settings.getStage3GraphFilename())) {
                cout << "Graph file " << Util::fileExists(settings.getStage3GraphFilename()) <<
                " not found. Aborting Stage 4" << endl;
                return;
        }
        
        string ifn = settings.getStage3GraphFilename();
        cout << "Loading graph from file: " << ifn << "..."; cout.flush();
        Util::startChrono();
        dBG.loadBinary(ifn);
        cout << "\n\tLoaded graph with " << dBG.getNumValidNodes() << " nodes and "
        << dBG.getNumValidArcs() << " arcs (" << Util::stopChronoStr() << ")\n";
        
        
        // Load or compute de node and edge models
        vector<double> wN(settings.getNumComp(), 1.0);
        vector<double> wE(settings.getNumComp(), 1.0);
        
        double initErrCov = settings.getInitErrCov();
        double initCov = settings.getInitCov();
        double initODF = settings.getInitODF();
        
        if (initCov < 0.0)
                initCov = dBG.getInitialKmerCovEstimate(initErrCov, 0.01);
        
        
        
        CovModel edgeModel(initErrCov, initODF, initCov, initODF, wE);
        CovModel nodeModel(initErrCov, initODF, initCov, initODF, wN);
        
        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMargin(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());
        
        NodeMap<Multiplicity> nodeMult(dBG.getNumValidNodes());
        EdgeMap<Multiplicity> edgeMult(dBG.getNumValidArcs());
        vector<NodeRep> nodes = dBG.getNodeReps(dBG.getNumValidNodes());
        vector<EdgeRep> edges = dBG.getEdgeReps(dBG.getNumValidArcs());
        
        // Get/train model
        if (Util::fileExists(settings.getNodeModelFilename(true)) && Util::fileExists(settings.getEdgeModelFilename(true))) {
                nodeModel = CovModel(settings.getNodeModelFilename(true));
                edgeModel = CovModel(settings.getEdgeModelFilename(true));
        } else {
                cout << "\nStage 3 model files not found, retraining model...\nInitial coverage estimate: " << fixed << initCov << endl;
                
                cout.precision(2);
                
                trainModel(dBG, settings,
                           nodeModel, edgeModel, true);
        }
        
        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;
        
        // Get all multiplicities
        cout << "Computing all multiplicities using " 
        << ((settings.approxInf()) ? "approximate" : ("exact (nb-size " + to_string(settings.getCRFDepth()) + ")"))
        << " inference" << endl;
        
        computeMultiplicities(dBG, settings,
                              nodeMult, edgeMult,
                              nodeModel, edgeModel,
                              settings.approxInf(), settings.computeMAP());
        
        cout << "Multiplicity determination finished in " << Util::stopChronoStr() << "\nGenerating output files... " << endl;
                
        if (Util::fileExists("genome.fasta")) {
                KmerNPPTable table(dBG, true);
                dBG.createKmerNPPTable(table);
                
                RefComp refComp(dBG, settings, table, "genome.fasta");
                NodeMap<int> trueNodeMult; EdgeMap<int> trueEdgeMult;
                refComp.getTrueMultiplicity(trueNodeMult, trueEdgeMult);
                
                ofstream nodeOFS(settings.getTrueNodeMultFilename());
                ofstream edgeOFS(settings.getTrueEdgeMultFilename());
                
                for (const auto& it : trueNodeMult)
                        nodeOFS << it.first << "\t" << it.second << "\n";
                
                // write the true edge multiplicities
                for (const auto& it : trueEdgeMult)
                        edgeOFS << it.first.getSrcID() << "\t"
                        << it.first.getDstID() << "\t"
                        << it.second << "\n";
                
                nodeOFS.close();
                edgeOFS.close();
        }
        
        // output files (estimated multiplicities for nodes/edges)
        ofstream nodeOFS(settings.getEstNodeMultFilename());
        ofstream edgeOFS(settings.getEstEdgeMultFilename());
        
        for (size_t i = 0; i < nodes.size(); i++) {
                SSNode node = dBG.getSSNode(nodes[i]);
                nodeOFS << nodes[i].getNodeID() << "\t"
                << nodeMult[nodes[i]].getExpMult() << "\t"
                << nodeMult[nodes[i]].getExpMultLProb() << "\t"
                << node.getAvgCov() << "\t"
                << node.getMarginalLength() << "\n";
        }
        
        for (size_t i = 0; i < edges.size(); i++) {
                edgeOFS << edges[i].getSrcID() << "\t"
                << edges[i].getDstID() << "\t"
                << edgeMult[edges[i]].getExpMult() << "\t"
                << edgeMult[edges[i]].getExpMultLProb() << "\t"
                << dBG.getArc(edges[i]).getCov() << "\n";
        }
                
        // Output a cytoscape graph for visual inspection
        drawFullCytoGraph(settings, dBG);
        
        cout << "... finished." << endl;
}


int main(int argc, char** argv)
{
        try {
                Settings settings(argc, argv);
                LibraryContainer libraries(settings.getReadFilename());

                cout << "Welcome to Detox version " << MYTOX_MAJOR_VERSION << "." << MYTOX_MINOR_VERSION << "." << MYTOX_PATCH_LEVEL;

#ifdef DEBUG
                cout << " (debug mode)" << endl;
#else
                cout << " (release mode)" << endl;
#endif
		// Read DBG and compute abundancies from read data
                stageOne(settings, libraries);
		// Compute coverage distributions for different repeat levels
                stageTwo(settings);

		/*
                if(settings.getVisGraphNode() > 0) {
                        visualiseSubgraph(settings);
                } else {
                        stageThreeCorrectGraph(settings);
                        StageFourComputeMult(settings);
                }
		*/

		stageAssemble(settings);
		
        } catch (exception& e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }

        cout << "Exiting... bye!" << endl;
        return EXIT_SUCCESS;
}
