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

#include <cstdlib>
#include <numeric>

#include "dbgraph.h"
#include "settings.h"
#include "refcomp.h"
#include "crfmult.h"
#include "readaln.h"
#include "threshold.h"

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


void stageThreshold(Settings& settings)
{
  cout << "\nEntering stage thresholding\n";
  cout << "===========================\n" << endl;
  
  if (Util::fileExists(settings.getThresholdFilename())) {
    cout << "File " << settings.getThresholdFilename()
	 << " exists. Skipping threshold computation.\n";
    return;
  }

  Util::startChrono();
  Thresholds th;

  th.computeThresholds("histogram.node.st2.dat", settings.getMisassLh());

  th.writeThresholds(settings.getThresholdFilename(), settings.getProbFilename());
  
  cout << "Stage thresholding finished in " << Util::stopChronoStr() << endl;	
}

bool compareNodes(SSNode n1, SSNode n2) {
  return (n1.getAvgCov() < n2.getAvgCov());
}

void stageAssemble(Settings& settings)
{
  cout << "\nEntering stage assemble\n";
  cout << "=======================\n" << endl;

  if (!Util::fileExists(settings.getThresholdFilename())) {
    cout << "File " << settings.getThresholdFilename()
	 << " does not exist. Compute the thresholds file and rerun.\n";
    return;
  }

  if (!Util::fileExists(settings.getProbFilename())) {
    cout << "File " << settings.getProbFilename()
	 << " does not exist. Compute the thresholds file and rerun.\n";
    return;
  }

  
  // Read the thresholds
  Thresholds th;
  th.readThresholds(settings.getThresholdFilename(), settings.getProbFilename());
  
  DBGraph dBG(settings);
  string ifn = settings.getStage1GraphFilename();
  cout << "Loading graph from file: " << ifn << "..."; cout.flush();
  Util::startChrono();
  dBG.loadBinary(ifn);
  
  cout << "\n\tLoaded " << dBG.getNumNodes() << " nodes and "
       << dBG.getNumArcs() << " arcs (" << Util::stopChronoStr() << ")\n";

  Util::startChrono();
  
  double nodeCutoff= 4.9;
  vector<NodeRep> nodeReps = dBG.getLowCovNodes(nodeCutoff);
  cout << "Selected " << nodeReps.size()
       << " nodes with coverage <= " << nodeCutoff << endl;

  dBG.removeNodes(nodeReps);
  dBG.concatenateNodes();

  // Examine each arc and remove it if it does not meet the threshold
  std::vector<EdgeRep> edges = dBG.getEdgeReps(dBG.getNumValidArcs());
  vector<EdgeRep> edgesToRemove;
  for (size_t i = 0; i < edges.size(); i++) {
    NodeID srcID = edges[i].getSrcID();
    NodeID dstID = edges[i].getDstID();
    Arc arc = dBG.getArc(edges[i]);
    if (arc.isValid()) {
      SSNode n = dBG.getSSNode(srcID);
#ifndef AVG_COV
      if (arc.getCov() >= th.get((int)n.getCount(n.getMarginalLength()-1))) {
	SSNode n2 = dBG.getSSNode(dstID);
	if (arc.getCov() >= th.get((int)n2.getCount(n2.getMarginalLength()-1))) {
	} else {
	  edgesToRemove.push_back(edges[i]);
	}
      } else {
	// These are to be removed
	edgesToRemove.push_back(edges[i]);
	/*
	  std::cout << arc.getCov() << std::endl;
	  std::cout << " " << dBG.getSSNode(srcID).getAvgCov() << std::endl;
	  std::cout << " " << dBG.getSSNode(dstID).getAvgCov() << std::endl;
	*/
      }
#else
      if (arc.getCov() >= th.get((int)n.getAvgCov())) {
      } else {
	// These are to be removed
	edgesToRemove.push_back(edges[i]);
	/*
	  std::cout << arc.getCov() << std::endl;
	  std::cout << " " << dBG.getSSNode(srcID).getAvgCov() << std::endl;
	  std::cout << " " << dBG.getSSNode(dstID).getAvgCov() << std::endl;
	*/
      }
#endif
    }
  }
  dBG.removeEdges(edgesToRemove);
  dBG.concatenateNodes();

  cout << "\tGraph has " << dBG.getNumValidNodes()
       << " nodes and " <<  dBG.getNumValidArcs() << " arcs\n";
  /*
  std::vector<EdgeRep> edges2 = dBG.getEdgeReps(dBG.getNumValidArcs());
  for (size_t i = 0; i < edges2.size(); i++) {
    NodeID srcID = edges2[i].getSrcID();
    NodeID dstID = edges2[i].getDstID();
    Arc arc = dBG.getArc(edges2[i]);
    if (arc.isValid()) {
      std::cout << arc.getCov() << " " << th.get((int)dBG.getSSNode(srcID).getAvgCov()) <<  std::endl;
      std::cout << " " << dBG.getSSNode(srcID).getAvgCov() << std::endl;
      std::cout << " " << dBG.getSSNode(dstID).getAvgCov() << std::endl;
    }
  }
  */

  // Generate unitigs
  dBG.defrag();
  dBG.setAllFlags1(false);
  vector<SSNode> nodes;
  nodes.resize(dBG.getNumNodes());
  for(int i = 1; i <= dBG.getNumNodes(); i++) {
    nodes[i-1] = dBG.getSSNode(i);
  }
  sort(nodes.begin(), nodes.end(), compareNodes);

  int id = 1;
  ofstream ofs("output.fa");
  ofstream ofsp("output.prob");
  int breaks = 0;
  for(int i = 0; i < dBG.getNumNodes(); i++) {
    SSNode start = nodes[i];
    if (!start.isValid())
      continue;
    //std::cout << start.getAvgCov() << std::endl;
    if (!start.getFlag1() && start.getAvgCov() >= th.get((int)start.getAvgCov())) {
      start.setFlag1(true);
      vector<NodeID> nodeSeq;
      string contig;
      vector<double> fwdP;
      vector<double> bwdP;
      nodeSeq.clear();
      nodeSeq.push_back(start.getNodeID());
      
      // Extend right
      SSNode current = start;
      while(current.numRightArcs() == 1) {
	NodeID rightID = current.rightBegin()->getNodeID();
	SSNode right = dBG.getSSNode(rightID);
	// don't merge palindromic repeats / loops
	if (right.getFlag1())
	  break;
	nodeSeq.push_back(rightID);
	right.setFlag1(true);
	current = right;
      }

      // Extend left
      std::reverse(nodeSeq.begin(), nodeSeq.end());
      for(int i = 0; i < nodeSeq.size(); i++) {
	nodeSeq[i] = -nodeSeq[i];
      }
      current = dBG.getSSNode(nodeSeq[nodeSeq.size()-1]);
      while(current.numRightArcs() == 1) {
	NodeID rightID = current.rightBegin()->getNodeID();
	SSNode right = dBG.getSSNode(rightID);
	// don't merge palindromic repeats / loops
	if (right.getFlag1())
	  break;
	nodeSeq.push_back(rightID);
	right.setFlag1(true);
	current = right;
      }

      /*
      for(auto it = nodeSeq.begin(); it != nodeSeq.end(); ++it) {
	std::cout << " " << *it << " (" << dBG.getSSNode(*it).getAvgCov() << ")";
      }
      std::cout << endl;
      */

#ifndef AVG_COV
      // Write out the contig(s)
      contig.clear();
      fwdP.clear();
      bwdP.clear();
      for(int i = 0; i < Kmer::getK()-1; i++) {
	fwdP.push_back(0.0);
      }
      for(size_t i = 0; i < nodeSeq.size(); i++) {
	size_t start = 0;
	SSNode n = dBG.getSSNode(nodeSeq[i]);
	if (i == 0) {
	  contig = n.getSequence().substr(start, Kmer::getK()-1);
	}
	for(size_t end = 0; end < n.getMarginalLength()-1; end++) {
	  if (n.getArcCount(end) < th.get(n.getCount(end)) || n.getArcCount(end) < th.get(n.getCount(end+1))) {
	    breaks++;
	    contig.append(n.getSequence().substr(start+Kmer::getK()-1, end-start+1));
	    if (contig.length() > Kmer::getK()) {
	      for(int j = 0; j < Kmer::getK()-1; j++) {
		bwdP.push_back(0.0);
	      }
	      if (contig.length()-1 != fwdP.size()) {
		std::cout << "Mismatch in contig length and fwdP size: " <<  contig.length() << " " << fwdP.size() << " " << "contig_" << id << std::endl; 
	      }
	      if (contig.length()-1 != bwdP.size()) {
		std::cout << "Mismatch in contig length and bwdP size: " <<  contig.length() << " " << bwdP.size() << " " << "contig_" << id << std::endl; 
	      }
	      ofs << ">contig_" << id << "\n";
	      Util::writeSeqWrap(ofs, contig, 60);
	      ofsp << ">contig_" << id << "_fwd\n";
	      Util::writeProbWrap(ofsp, fwdP, 60);
	      ofsp << ">contig_" << id << "_bwd\n";
	      Util::writeProbWrap(ofsp, bwdP, 60);
	      id++;
	    }
	    start = end+1;
	    contig.clear();
	    fwdP.clear();
	    bwdP.clear();
	    for(int j = 0; j < Kmer::getK()-1; j++) {
	      fwdP.push_back(0.0);
	    }
	    contig = n.getSequence().substr(start, Kmer::getK()-1);
	  } else {
	    fwdP.push_back(th.getProb(n.getCount(end), n.getArcCount(end)));
	    //std::cout << "internal: " << n.getCount(end) << " " << n.getArcCount(end) << " " << th.getProb(n.getCount(end), n.getArcCount(end)) << std::endl;
	    bwdP.push_back(th.getProb(n.getCount(end+1), n.getArcCount(end)));
	  }
	}
	contig.append(n.getSequence().substr(start+Kmer::getK()-1));
	if (i < nodeSeq.size()-1){
	  fwdP.push_back(th.getProb(n.getCount(n.getMarginalLength()-1), n.rightArc(nodeSeq[i+1])->getCov()));
	  //std::cout << "external: " << n.getCount(n.getMarginalLength()-1) << " " << n.rightArc(nodeSeq[i+1])->getCov() << " " << th.getProb(n.getCount(n.getMarginalLength()-1), n.rightArc(nodeSeq[i+1])->getCov()) << std::endl;
	  SSNode next = dBG.getSSNode(nodeSeq[i+1]);
	  bwdP.push_back(th.getProb(next.getCount(0), n.rightArc(nodeSeq[i+1])->getCov()));
	}
      }
      if (contig.length() > Kmer::getK()) {
	for(int j = 0; j < Kmer::getK()-1; j++) {
	  bwdP.push_back(0.0);
	}
	if (contig.length()-1 != fwdP.size()) {
	  std::cout << "Mismatch in contig length and fwdP size: " <<  contig.length() << " " << fwdP.size() << " " << "contig_" << id << std::endl; 
	}
	if (contig.length()-1 != bwdP.size()) {
	  std::cout << "Mismatch in contig length and bwdP size: " <<  contig.length() << " " << bwdP.size() << " " << "contig_" << id << std::endl; 
	}
	ofs << ">contig_" << id << "\n";
	Util::writeSeqWrap(ofs, contig, 60);
	ofsp << ">contig_" << id << "_fwd\n";
	Util::writeProbWrap(ofsp, fwdP, 60);
	ofsp << ">contig_" << id << "_bwd\n";
	Util::writeProbWrap(ofsp, bwdP, 60);
	id++;
      }
	
      
#else
      // Write out the contig
      dBG.convertNodesToString(nodeSeq, contig);
      ofs << ">contig_" << id << "\n";
      Util::writeSeqWrap(ofs, contig, 60);      
      id++;
#endif
    }
  }
  std::cout << "Number of breaks due to coverage: " << breaks << std::endl;
  ofs.close();
  ofsp.close();
  
  cout << "Stage assemble finished in " << Util::stopChronoStr() << endl;

  //  dBG.writeContigs("output2.fa");

}


int main(int argc, char** argv)
{
        try {
                Settings settings(argc, argv);
                LibraryContainer libraries(settings.getReadFilename());

                cout << "Welcome to SAMA version " << SAMA_MAJOR_VERSION << "." << SAMA_MINOR_VERSION << "." << SAMA_PATCH_LEVEL;

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

		stageThreshold(settings);
		
		stageAssemble(settings);
		
        } catch (exception& e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }

        cout << "Exiting... bye!" << endl;
        return EXIT_SUCCESS;
}
