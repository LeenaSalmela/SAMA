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

#ifndef CRFMULT_H
#define CRFMULT_H

#include <vector>
#include <cstdlib>
#include <map>
#include <queue>
#include <atomic>

#include "coverage.h"
#include "pgm/factor.h"
#include "dbgraph.h"
#include "bitvec.h"
#include "libdai/dai/factorgraph.h"
#include "libdai/dai/factor.h"

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class DBGraph;
class WorkLoadBalancer;
class NodeRep;
class EdgeRep;

// ============================================================================
// DIJKSTRA AUXILIARY CLASSES
// ============================================================================

class PathDFS {

public:
        NodeID nodeID;          // current node identifier
        NodeID prevID;          // previous node identifier
        size_t length;          // length to current node

        /**
         * Default constructor
         * @param nodeID node identifier
         * @param prevID previous node identifier
         * @param length length to the current node
         */
        PathDFS(NodeID nodeID, NodeID prevID, size_t length) :
                nodeID(nodeID), prevID(prevID), length(length) {};
};

struct PathDFSComp {

        /**
         * Compare two paths (by length) (strict weak ordering (irreflexive))
         */
        bool operator()(const PathDFS& f, const PathDFS& s) {
                return f.length > s.length;
        }
};

// ============================================================================
// CRF PERFORMANCE COUNTER
// ============================================================================

class CRFPerfCounter {

public:
        size_t totNumCRF;
        size_t totNumNodes;
        size_t totNumEdges;

        /**
         * Initialization constructor
         */
        CRFPerfCounter() : totNumCRF(0), totNumNodes(0), totNumEdges(0) {}

        /**
         * Merge counts
         * @param rhs Right hand size
         **/
        void merge(const CRFPerfCounter& rhs) {
                totNumCRF += rhs.totNumCRF;
                totNumNodes += rhs.totNumNodes;
                totNumEdges += rhs.totNumEdges;
        }

        /**
         * Get the average number of nodes in a CRF neighborhood
         * @return The average number of edges in a CRF neighborhood
         **/
        double getAvgNumNodes() const {
                return (double)totNumNodes / totNumCRF;
        }

        /**
         * Get the average number of edges in a CRF neighborhood
         * @return The average number of edges in a CRF neighborhood
         **/
        double getAvgNumEdges() const {
                return (double)totNumEdges / totNumCRF;
        }
};

// ============================================================================
// CONDITIONAL RANDOM FIELDS SOLVER (PER THREAD)
// ============================================================================

/**
 * This class is designed for speed. Most of the containers are allocated as
 * class members instead of local variables to avoid repeated (de)allocation
 * of these containers over many multiplicity computations.
 * Therefore, some functions may have side effects: they modify class members
 * containers. This class is therefore not thread-safe. The idea is to run
 * a single class object within each thread.
 */

class CRFSolver {

private:
        const DBGraph& dBG;             // const-ref to de Bruijn graph

        int multMargin;                 // number of alt mult (one-sided)
        size_t maxFactorSize;           // maximum size of intermediate factor
        double flowStrength;            // value in flow-cpd when mults do not agree

        // containers below are class members to allow their reuse over many
        // multiplicity computations, thus avoiding repeated (de)allocation
        std::vector<NodeRep> nodes;     // nodes in subgraph
        std::vector<EdgeRep> edges;     // edges in subgraph
        std::priority_queue<NodeDFS, std::vector<NodeDFS>, NodeDFSComp> pq;
        Bitvec bitvec;

        /**
         * Get a singleton factor with a multinomial over multiplicities
         * @param varID Variable identifier
         * @param card Variable cardinality
         * @param firstMult Multiplicity of the first value
         * @param covModel Coverage model
         * @param numObs Number of observations
         * @param numIndepOb FIXME !!!
         * @return Singleton multiplicity factor
         */
        static Factor createSingletonFactor(int varID, int card, int firstMult,
                                            const CovModel& covModel,
                                            double numObs, int numIndepObs);

        /**
         * Get a Libdai format singleton factor with a multinomial over multiplicities
         * @param varID Variable identifier
         * @param varcard Variable cardinality
         * @param firstMult Multiplicity of the first value
         * @param covModel Coverage model
         * @param numObs Number of observations
         * @param numIndepOb FIXME !!!
         * @return Singleton multiplicity factor
         */
        static dai::Factor createLDSingletonFactor(int varID, int varCard, int firstMult,
                                                   const CovModel& covModel,
                                                   const double numObs, int numIndepObs);


        /**
         * Get a flow-conservation factor
         * @param sumVarID Variable identifier of the sum
         * @param sumCard Cardinality of the sum variable
         * @param sumFirstMult Multiplicity of the first sum value
         * @param termVarID Variable identifiers of the terms
         * @param termCard Cardinality of the term variables
         * @param termFirstMult Multiplicity of the first term values
         * @param palindromic Is the edge palindromic
         * @return Flow-conservation factor
         */
        static Factor createFlowFactor(int sumVarID, int sumCard, int sumFirstMult,
                                       const std::vector<int>& termVarID,
                                       const std::vector<int>& termCard,
                                       const std::vector<int>& termFirstMult,
                                       const std::vector<bool>& palindromic,
                                       double flowStrength);
        
        /**
         * Get a Libdai format flow-conservation factor
         * @param sumVarID Variable identifier of the sum
         * @param sumCard Cardinality of the sum variable
         * @param sumFirstMult Multiplicity of the first sum value
         * @param termVarID Variable identifiers of the terms
         * @param termCard Cardinality of the term variables
         * @param termFirstMult Multiplicity of the first term values
         * @param palindromic Is the edge palindromic
         * @return Flow-conservation factor
         */
        static dai::Factor createLDFlowFactor(int sumVarID, int sumCard, int sumFirstMult,
                                       const std::vector<int>& termVarID,
                                       const std::vector<int>& termCard,
                                       const std::vector<int>& termFirstMult,
                                       const std::vector<bool>& palindromic,
                                       double flowStrength);
        

        /**
         * Internal function to find a CRF-compliant subgraph.
         * Output is written to nodes and edges as a side-effect.
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(int maxDepth);

        /**
         * Get a subgraph of the de Bruijn graph centered around a seed node
         * Output is written to nodes and edges as a side-effect
         * @param seedNode Seed node representative
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(NodeRep seedNode, int maxDepth);

        /**
         * Get a subgraph of the de Bruijn graph centered around a seed edge
         * Output is written to nodes and edges as a side-effect
         * @param seedNode Seed node representative
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(EdgeRep seedEdge, int maxDepth);

        /**
         * Compute the multiplicity using a CRF model
         * @param node2var Mapping between dBG nodes and CRF variables
         * @param edge2var Mapping between dBG edges and CRF variables
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param resMult Resulting multiplicity (output)
         * @param targetCRFVar CRF variable to keep
         */
        bool solveSubgraph(const std::map<NodeRep, int>& node2var,
                           const std::map<EdgeRep, int>& edge2var,
                           const CovModel& nodeCovModel,
                           const CovModel& edgeCovModel,
                           Multiplicity& resMult, int targetCRFVar) const;

public:
        CRFPerfCounter perfCounter;     // performance counter

        /**
         * Constructor
         * @param dBG Const-reference to de Bruijn graph
         * @param multMargin Number of alternative multiplicities (one-sided)
         * @param maxFactorSize Maximum size of intermediate factor
         * @param numThreads Number of threads
         * @param threadWork Amount of work (nodes/edges) per work chunk
         */
        CRFSolver(const DBGraph& dBG, int multMargin, size_t maxFactorSize,
                  double flowStrength) : dBG(dBG), multMargin(multMargin),
                  maxFactorSize(maxFactorSize), flowStrength(flowStrength),
                  bitvec(dBG.getNumNodes()+1) {}
                  
        CRFSolver(const DBGraph& dBG, int multMargin, size_t maxFactorSize,
                  double flowStrength, const std::vector<NodeRep>& n, const std::vector<EdgeRep>& e) : dBG(dBG), multMargin(multMargin),
                            maxFactorSize(maxFactorSize), flowStrength(flowStrength), nodes(n), edges(e),
                            bitvec(dBG.getNumNodes()+1) {}

        /**
         * Compute the node multiplicity using a CRF model
         * @param node Node
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param graphDepth CRF graph depth (in terms of number of nodes)
         * @return true of false (when flow is not conserved)
         */
        bool checkFlow(NodeRep node, const CovModel& nodeCovModel,
                       const CovModel& edgeCovModel, int graphDepth);

        /**
         * Compute the node multiplicity using a CRF model
         * @param node Node
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param graphDepth CRF graph depth (in terms of number of nodes)
         */
        Multiplicity getNodeMultiplicity(NodeRep node,
                                         const CovModel& nodeCovModel,
                                         const CovModel& edgeCovModel,
                                         int graphDepth);

        /**
         * Compute the edge multiplicity using a CRF model
         * @param edge Edge
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param graphDepth CRF graph depth (in terms of number of nodes)
         */
        Multiplicity getEdgeMultiplicity(EdgeRep edge,
                                         const CovModel& nodeCovModel,
                                         const CovModel& edgeCovModel,
                                         int graphDepth);
        
        /**
         * Create a factor graph for the full dBG
         * 
         */
        dai::FactorGraph getLibdaiFG(const std::map<NodeRep, int>& node2var,
                                     const std::map<EdgeRep,int>& edge2var,
                                     std::vector<int>& firstMult,
                                     const CovModel& nodeCovModel,
                                     const CovModel& edgeCovModel);
        
	/*
	 * Approximate inference on a subgraph to be used for the cytoscape graph
	 * Pass along the nodes and edges ID vector to determine correct order for node and edge Mult
	 * because we need the outer edges for the crf, the subgraph will be extracted anew: 
	 * thus central node and required graph depth need to be passed again.
	 */
	void approxSubgraphMult(NodeRep node,
				NodeMap<Multiplicity>& nodeMult,
				EdgeMap<Multiplicity>& edgeMult,
				const CovModel& nodeCovModel,
				const CovModel& edgeCovModel,
				int graphDepth, bool MAP);
        /**
         * Get subgraph from the de Bruijn graph 
         * that contains all nodes reachable from seedNode
         * @param seedNode Seed node representative
         */
        std::vector<NodeID> getFullConnectedComponent(NodeRep seedNode);
        
        /**
         * Use approximate inference to compute te multiplicities of the passed nodes
         * and arcs. Assumes that these form a connected subgraph
         */
        void approxMult(NodeMap<Multiplicity>& nodeMult,
                        const CovModel& nodeCovModel,
                        EdgeMap<Multiplicity>& edgeMult,
                        const CovModel& edgeCovModel, 
                        std::string modelID,
                        bool MAP, bool visualise = false);
};

// ============================================================================
// ESTIMATE NODE/EDGE MULTIPLICITIES USING CONDITIONAL RANDOM FIELDS
// ============================================================================

class CRFMult {

private:
        const DBGraph& dBG;             // const-ref to de Bruijn graph

        int maxGraphDepth;              // maximum depth of the subgraph
        int multMargin;                 // number of alt mult (one-sided)
        size_t maxFactorSize;           // maximum size of intermediate factor
        double flowStrength;            // value in flow-cpd when mults do not agree

        size_t numThreads;              // number of threads
        size_t threadWork;              // amount of work per thread

        mutable CRFPerfCounter totPerfCounter;  // total performance counter

        /**
         * Compute the node conservation of flow
         * @param solver Per-thread solver
         * @param wlb Reference to the workload balancer
         * @param nodes Nodes to handle
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param flowOK Flag to indicate whether or not the flow is OK
         */
        void checkNodeFlow(CRFSolver& solver,
                           WorkLoadBalancer& wlb,
                           const std::vector<NodeRep>& nodes,
                           const CovModel& nodeCovModel,
                           const CovModel& edgeCovModel,
                           std::vector<bool>& flowOK) const;

        /**
         * Compute the node multiplicities (thread entry)
         * @param wlb Reference to the workload balancer
         * @param nodes Nodes to handle
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param nodeMult Estimated node multiplicity (output)
         */
        void computeNodeMult(CRFSolver& solver,
                             WorkLoadBalancer& wlb,
                             const std::vector<NodeRep>& nodes,
                             const CovModel& nodeCovModel,
                             const CovModel& edgeCovModel,
                             NodeMap<Multiplicity>& nodeMult) const;

        /**
         * Compute the edge multiplicities (thread entry)
         * @param wlb Reference to the workload balancer
         * @param edges Edges to handle
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param edgeMult Estimated edge multiplicity (output)
         */
        void computeEdgeMult(CRFSolver& solver,
                             WorkLoadBalancer& wlb,
                             const std::vector<EdgeRep>& edges,
                             const CovModel& nodeCovModel,
                             const CovModel& edgeCovModel,
                             EdgeMap<Multiplicity>& edgeMult) const;

        /**
         * Compute the node/edge multiplicities and check for changes
         * @param nodeMult Estimated node multiplicity (input/output)
         * @param nodeCovModel Node coverage model
         * @param edgeMult Estimated edge multiplicity (input/output)
         * @param edgeCovModel Edge coverage model
         * @param epsilon Tolerance on fraction of altered nodes/edges
         * @return True if a fraction of > epsilon nodes or edges changed
         */
        bool Estep(NodeMap<Multiplicity>& nodeMult,
                   const CovModel& nodeCovModel,
                   EdgeMap<Multiplicity>& edgeMult,
                   const CovModel& edgeCovModel,
                   double epsilon, bool approxInf = false,
                   bool map = false) const;

        /**
         * Update the node/edge models given the multiplicities
         * @param nodeMult Estimated node multiplicity
         * @param nodeCovModel Node coverage model (input/output)
         * @param edgeMult Estimated edge multiplicity
         * @param edgeCovModel Edge coverage model (input/output)
         * @param epsilon Tolerance on the lambda_1 estimate
         * @return True if the lambda_1 component change is > epsilon
         */
        bool Mstep(const NodeMap<Multiplicity>& nodeMult,
                   CovModel& nodeCovModel,
                   const EdgeMap<Multiplicity>& edgeMult,
                   CovModel& edgeCovModel, double epsilon, 
                   bool fixedZero = false) const;

public:
        /**
         * Constructor
         * @param dBG Const-reference to de Bruijn graph
         * @param maxGraphDepth Maximum depth of a subgraph
         * @param multMargin Number of alternative multiplicities (one-sided)
         * @param maxFactorSize Maximum size of intermediate factor
         * @param numThreads Number of threads
         * @param threadWork Amount of work (nodes/edges) per work chunk
         */
        CRFMult(const DBGraph& dBG, int maxGraphDepth,
                int multMargin, size_t maxFactorSize, double flowStrength,
                size_t numThreads = 1, size_t threadWork = 1000) : dBG(dBG), maxGraphDepth(maxGraphDepth),
                multMargin(multMargin), maxFactorSize(maxFactorSize), flowStrength(flowStrength),
                numThreads(numThreads), threadWork(threadWork) {}

        /**
         * Check whether the flow around a set of nodes is OK
         * @param nodes Nodes to handle
         * @param flowOK Flag to indicate whether or not the flow is OK
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Edge coverage model
         */
        void checkFlow(const std::vector<NodeRep>& nodes,
                       std::vector<bool>& flowOK,
                       const CovModel& nodeCovModel,
                       const CovModel& edgeCovModel) const;

        /**
         * Compute the normalized node/edge multiplicities
         * @param nodeMult Estimated node multiplicity (input/output)
         * @param nodeCovModel Node coverage model
         * @param edgeMult Estimated edge multiplicity (input/output)
         * @param edgeCovModel Edge coverage model
         */
        void computeMult(NodeMap<Multiplicity>& nodeMult,
                         const CovModel& nodeCovModel,
                         EdgeMap<Multiplicity>& edgeMult,
                         const CovModel& edgeCovModel) const;

        /**
         * Uses approximate inference to compute all multiplicities in the dBG 
         * After computations, nodeMult and edgeMult will have entries for all valid nodes/arcs
         */
        void approxMultAll(NodeMap<Multiplicity>& nodeMult,
                           const CovModel& nodeCovModel,
                           EdgeMap<Multiplicity>& edgeMult,
                           const CovModel& edgeCovModel,
                           bool MAP=false, bool singleCRF=false) const;
                           
        /**
         * Compute the multiplicities and models using expectation-maximization
         * @param nodeMult Estimated normalized node multiplicity (input/output)
         * @param nodeCovModel Node coverage model (input/output)
         * @param edgeMult Estimated normalized edge multiplicity (input/output)
         * @param edgeCovModel Edge coverage model (input/output)
         * @param epsilon Tolerance on the lambda_1 estimate
         * @param maxIter Maximum number of iterations
         */
        int computeMultEM(NodeMap<Multiplicity>& nodeMult,
                          CovModel& nodeCovModel,
                          EdgeMap<Multiplicity>& edgeMult,
                          CovModel& edgeCovModel,
                          double epsilon, int maxIter,
                          bool approxInf = false, bool map = false,
                          bool fixedZero = false);

        /**
         * Get the CRF performance counter
         * @return The CRF performance counter
         */
        CRFPerfCounter getPerfCounter() const {
                return totPerfCounter;
        }
};

#endif
