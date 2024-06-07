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

#include <numeric>
#include <thread>

#include "crfmult.h"
#include "dbgraph.h"
#include "pgm/pgminference.h"
#include "libdai/dai/alldai.h"  // Include main libDAI header file
#include "libdai/dai/jtree.h"
#include "libdai/dai/bp.h"
#include "libdai/dai/parallelbp.h"
#include "libdai/dai/factor.h"
#include "libdai/dai/factorgraph.h"
#include "libdai/dai/var.h"
#include "libdai/dai/varset.h"
#include "arc.h"
#include "ssnode.h"
#include "arc.h"

using namespace std;

// ============================================================================
// CONDITIONAL RANDOM FIELDS SOLVER (PER THREAD)
// ============================================================================

void CRFSolver::getSubgraph(int maxDepth)
{
        while (!pq.empty()) {
                // get and erase the current node
                NodeDFS currTop = pq.top();
                pq.pop();
                NodeID thisID = currTop.nodeID;
                int thisDepth = currTop.depth;

                // if the node was already handled, skip
                if (bitvec[abs(thisID)])
                        continue;

                nodes.push_back(NodeRep(thisID));
                SSNode n = dBG.getSSNode(thisID);

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        if (bitvec[abs(rightID)])       // edge already added?
                                continue;

                        edges.push_back(EdgeRep(thisID, rightID));
                        if (thisDepth < maxDepth)
                                pq.push(NodeDFS(rightID, thisDepth+1));
                }

                // process the left arcs
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        if (bitvec[abs(leftID)])        // edge already added?
                                continue;
                        if (leftID == thisID)   // edge already added as right arc
                                continue;

                        edges.push_back(EdgeRep(leftID, thisID));
                        if (thisDepth < maxDepth)
                                pq.push(NodeDFS(leftID, thisDepth + 1));
                }

                // mark this node as handled
                bitvec[abs(thisID)] = true;
        }

        // reset all flags to false
        for (auto it : nodes)
                bitvec[it.getNodeID()] = false;
}

void CRFSolver::getSubgraph(NodeRep seedNode, int maxDepth)
{
        nodes.clear(); edges.clear();

        if (maxDepth == 0) {
                nodes.push_back(seedNode);
                return;
        }

        // add the seed node to the priority queue
        pq.push(NodeDFS(seedNode, 0));

        getSubgraph(maxDepth);
}

void CRFSolver::getSubgraph(EdgeRep seedEdge, int maxDepth)
{
        nodes.clear(); edges.clear();

        // special case (maxDepth == 0)
        if (maxDepth == 0) {
                edges.push_back(seedEdge);
                return;
        }

        // other cases (maxDepth > 0)
        pq.push(NodeDFS(seedEdge.getSrcID(), 1));
        pq.push(NodeDFS(seedEdge.getDstID(), 1));

        getSubgraph(maxDepth);
}

vector<NodeID> CRFSolver::getFullConnectedComponent(NodeRep seedNode)
{
        nodes.clear(); edges.clear();
        vector<NodeID> usedIDs;
        usedIDs.clear();
        
        // add the seed node to the priority queue
        pq.push(NodeDFS(seedNode, 0));
        
        // Same functionallity as getSubgraph but no maxDepth check
        while (!pq.empty()) {
                // get and erase the current node
                NodeDFS currTop = pq.top();
                pq.pop();
                NodeID thisID = currTop.nodeID;
                int thisDepth = currTop.depth;
                
                // if the node was already handled, skip
                if (bitvec[abs(thisID)])
                        continue;
                
                nodes.push_back(NodeRep(thisID));
                usedIDs.push_back(thisID);
                SSNode n = dBG.getSSNode(thisID);
                
                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        if (bitvec[abs(rightID)])       // edge already added?
                                continue;
                        
                        edges.push_back(EdgeRep(thisID, rightID));
                        pq.push(NodeDFS(rightID, thisDepth+1));
                }
                
                // process the left arcs
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        if (bitvec[abs(leftID)])        // edge already added?
                                continue;
                        if (leftID == thisID)   // edge already added as right arc
                                continue;
                        
                        edges.push_back(EdgeRep(leftID, thisID));
                        pq.push(NodeDFS(leftID, thisDepth + 1));
                }
                
                // mark this node as handled
                bitvec[abs(thisID)] = true;
        }
        
        // reset all flags to false
        for (auto it : nodes)
                bitvec[it.getNodeID()] = false;
        
        return usedIDs;
}

Factor CRFSolver::createSingletonFactor(int varID, int varCard, int firstMult,
                                        const CovModel& covModel,
                                        const double numObs, int numIndepObs)
{
        vector<int> var(1, varID);
        vector<int> card(1, varCard);
        vector<double> val(varCard);

        for (int i = 0; i < varCard; i++) {
                int multiplicity = firstMult + i;
                val[i] = covModel.getLogProb(numObs, multiplicity);
        }

        /*if ((firstMult == 0) && (val.size() > 1))
                val[0] = max(val[0], log(10e-3) + val[1]);*/

        for (int i = 0; i < varCard; i++)
                val[i] *= numIndepObs;

        return Factor(move(var), move(card), move(val));
}

dai::Factor CRFSolver::createLDSingletonFactor(int varID, int varCard, int firstMult,
                                               const CovModel& covModel,
                                               const double numObs, int numIndepObs)
{
        dai::VarSet vars(dai::Var(varID,varCard));
        vector<double> val(varCard);

        double m = numeric_limits<double>::lowest();
        for (int i = 0; i < varCard; i++) {
                int multiplicity = firstMult + i;
                val[i] = numIndepObs * covModel.getLogProb(numObs, multiplicity);
                m = max(val[i], m);
        }
        
        for (int i  = 0; i< varCard; i++)
                val[i] = max(exp(val[i] - m), numeric_limits<double>::min());

        return dai::Factor(move(vars), move(val));
}


Factor CRFSolver::createFlowFactor(int sumVarID, int sumCard, int sumFirstMult,
                                   const vector<int>& termVarID,
                                   const vector<int>& termCard,
                                   const vector<int>& termFirstMult,
                                   const vector<bool>& palindromic,
                                   double flowStrength)
{
        assert(termVarID.size() == termFirstMult.size());

        const double sumOK = flowStrength;
        const double sumNOK = 1.0;

        // create the variable vector
        vector<int> var = { sumVarID };
        var.insert(var.end(), termVarID.begin(), termVarID.end());

        // create the cardinality vector
        vector<int> card = { sumCard };
        card.insert(card.end(), termCard.begin(), termCard.end());

        // create the value vector
        size_t numVal = accumulate(card.begin(), card.end(), 1ull, multiplies<size_t>());
        vector<double> val(numVal);

        Assignment assignment(card);
        for (size_t i = 0; i < val.size(); i++) {
                int sumTerms = 0;
                for (size_t j = 0; j < termFirstMult.size(); j++) {
                        // count palindromic arcs double
                        int c = (palindromic[j]) ? 2 : 1;
                        sumTerms += c * (termFirstMult[j] + assignment[j+1]);
                }

                bool flag = ((assignment[0] + sumFirstMult) == sumTerms);
                val[i] = (flag) ? log(sumOK) : log(sumNOK);
                assignment.increment();
        }

        return Factor(move(var), move(card), move(val));
}

dai::Factor CRFSolver::createLDFlowFactor(int sumVarID, int sumCard, int sumFirstMult,
                                   const vector<int>& termVarID,
                                   const vector<int>& termCard,
                                   const vector<int>& termFirstMult,
                                   const vector<bool>& palindromic,
                                   double flowStrength)
{
        assert(termVarID.size() == termFirstMult.size());

        const double sumOK = flowStrength;
        const double sumNOK = 1.0;
	
	// sort edges (and their cards, firstmult and palindromic indicator) based on varID (libdai sorts them internally, but the value vector is created based on the order we pass here). NodeID will always be the smallest so that's okay
	vector<size_t> idx(termVarID.size());
	iota(idx.begin(), idx.end(), 0);
	sort(idx.begin(), idx.end(),
	     [&termVarID](size_t i1, size_t i2) {return termVarID[i1] < termVarID[i2];});
	vector<int> termVarIDsort, termCardSort, termFirstMultSort;
	vector<bool> palindromicSort;
	
	for (size_t i: idx){
		termVarIDsort.push_back(termVarID[i]);
		termCardSort.push_back(termCard[i]);
		termFirstMultSort.push_back(termFirstMult[i]);
		palindromicSort.push_back(palindromic[i]);
	}

        // create the variable vector
        vector<int> var = { sumVarID };
        var.insert(var.end(), termVarIDsort.begin(), termVarIDsort.end());

        // create the cardinality vector
        vector<int> card = { sumCard };
        card.insert(card.end(), termCardSort.begin(), termCardSort.end());
        
        dai::VarSet vars(dai::Var(sumVarID,sumCard));
        for (int i=0; i< termVarID.size(); i++) {
                vars.insert(dai::Var(termVarIDsort[i], termCardSort[i]));
        }

        // create the value vector
        size_t numVal = accumulate(card.begin(), card.end(), 1ull, multiplies<size_t>());
        vector<double> val(numVal);

        Assignment assignment(card);
        for (size_t i = 0; i < val.size(); i++) {
                int sumTerms = 0;
                for (size_t j = 0; j < termFirstMultSort.size(); j++) {
                        // count palindromic arcs double
                        int c = (palindromicSort[j]) ? 2 : 1;
                        sumTerms += c * (termFirstMultSort[j] + assignment[j+1]);
                }

                bool flag = ((assignment[0] + sumFirstMult) == sumTerms);
                val[i] = (flag) ? sumOK : sumNOK;
                assignment.increment();
        }

        return dai::Factor(move(vars), move(val));
}


bool CRFSolver::solveSubgraph(const map<NodeRep, int>& node2var,
                              const map<EdgeRep, int>& edge2var,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              Multiplicity& resMult, int targetCRFVar) const
{
        // create the singleton factors
        list<Factor> factorList;
        vector<int> card(node2var.size() + edge2var.size());
        vector<int> firstMult(node2var.size() + edge2var.size());

        for (const auto it : node2var) {        // first the nodes
                // always compute firstMult and card even if factor is not added
                SSNode node = dBG.getSSNode(it.first.getNodeID());
                const int& nodeVarID = it.second;
                int expMult = nodeCovModel.getExpMult(node.getAvgCov());
                firstMult[nodeVarID] = max(0, expMult - multMargin);
                card[nodeVarID] = expMult + multMargin - firstMult[nodeVarID] + 1;

                // add a singleton node factor in case:
                // - the node is bigger than 2k (less correlation with edges) OR
                // - there are no edges (in case of isolated nodes) OR
                // - there is a single incoming edge OR
                // - there is a single outgoing edge
                if ((node.getMarginalLength() < 2*Kmer::getK()) &&
                    (node.numLeftArcs() > 1) && (node.numRightArcs() > 1) && ! (edge2var.size() == 0))
                        continue;

                int nIndptObs = max<int>(1, node.getMarginalLength() / (2*Kmer::getK()));
                Factor F = createSingletonFactor(nodeVarID, card[nodeVarID],
                                                 firstMult[nodeVarID],
                                                 nodeCovModel, node.getAvgCov(), nIndptObs);
                factorList.push_back(F);
                
                //cout << "Node: " << it.first.getNodeID() << "\nFirstMult: " << firstMult[nodeVarID] << "\n" << F << endl;

                /*if (abs(it.first.getNodeID()) == 1202) {
                        cout << "Node: " << F << endl;
                }*/
        }

        for (const auto it : edge2var) {        // then the edges
                // compute firstMult and card even if factor is not added
                const EdgeRep& edge = it.first;
                const int& edgeVarID = it.second;
                double cov = dBG.getArc(edge).getCov();
                int expMult = edgeCovModel.getExpMult(cov);
                firstMult[edgeVarID] = max(0, expMult - multMargin);
                card[edgeVarID] = expMult + multMargin - firstMult[edgeVarID] + 1;

                // add a singleton edge factor in case:
                // - an adjacent node is not in the CRF (extremal edge) OR
                // - both adjacent nodes have multiple edges on that side
                if ((node2var.find(NodeRep(edge.getSrcID())) != node2var.end()) &&
                    (node2var.find(NodeRep(edge.getDstID())) != node2var.end()) && (! (node2var.size() == 0))) {
                        SSNode src = dBG.getSSNode(edge.getSrcID());
                        if (src.numRightArcs() == 1)
                                continue;
                        SSNode dst = dBG.getSSNode(edge.getDstID());
                        if (dst.numLeftArcs() == 1)
                                continue;
                }

                Factor F = createSingletonFactor(edgeVarID, card[edgeVarID],
                                                 firstMult[edgeVarID],
                                                 edgeCovModel, cov, 1);
                factorList.push_back(F);
                //cout << "Edge: " << edge << "\n" << F << endl;

                /*if (edge == EdgeRep(1202, -3354)) {
                        cout << "Edge: " << F << endl;
                }*/
        }

        // create the flow conservation factors
        if ( !edge2var.empty() ){
        for (const auto it : node2var) {
                NodeID currID = it.first.getNodeID();
                const int& nodeVarID = it.second;
                SSNode node = dBG.getSSNode(currID);

                // left flow conservation factor
                vector<int> lEdgeVarID, lEdgeVarCard, lEdgeFirstMult;
                vector<bool> lEdgePalindromic;
                for (ArcIt lIt = node.leftBegin(); lIt != node.leftEnd(); lIt++) {
                        EdgeRep edge(lIt->getNodeID(), currID);
                        int edgeVarID = edge2var.at(edge);

                        lEdgeVarID.push_back(edgeVarID);
                        lEdgeVarCard.push_back(card[edgeVarID]);
                        lEdgeFirstMult.push_back(firstMult[edgeVarID]);
                        lEdgePalindromic.push_back(lIt->getNodeID() == -currID);
                }

                if (!lEdgeVarID.empty()) {
                        Factor F = createFlowFactor(nodeVarID, card[nodeVarID],
                                                    firstMult[nodeVarID], lEdgeVarID,
                                                    lEdgeVarCard, lEdgeFirstMult,
                                                    lEdgePalindromic, flowStrength);
                        factorList.push_back(F);
                }


                // right flow conservation factor
                vector<int> rEdgeVarID, rEdgeVarCard, rEdgeFirstMult;
                vector<bool> rEdgePalindromic;
                for (ArcIt rIt = node.rightBegin(); rIt != node.rightEnd(); rIt++) {
                        EdgeRep edge(currID, rIt->getNodeID());
                        int edgeVarID = edge2var.at(edge);

                        rEdgeVarID.push_back(edgeVarID);
                        rEdgeVarCard.push_back(card[edgeVarID]);
                        rEdgeFirstMult.push_back(firstMult[edgeVarID]);
                        rEdgePalindromic.push_back(currID == -rIt->getNodeID());
                }

                if (!rEdgeVarID.empty()) {
                        Factor F = createFlowFactor(nodeVarID, card[nodeVarID],
                                                    firstMult[nodeVarID], rEdgeVarID,
                                                    rEdgeVarCard, rEdgeFirstMult,
                                                    rEdgePalindromic, flowStrength);
                        factorList.push_back(F);
                }
        }
        }

        // Variable elimination
        bool retVal = PGMInference::solveVE(factorList, targetCRFVar, maxFactorSize);
        if (!retVal)
                return false;

        resMult = Multiplicity(firstMult[targetCRFVar],
                               factorList.front().getVal());
        resMult.normalize();

        //cout << "RESULT: " << resMult << endl;

        return true;
}

bool CRFSolver::checkFlow(NodeRep node, const CovModel& nodeCovModel,
                          const CovModel& edgeCovModel, int graphDepth)
{
        getSubgraph(node, graphDepth);

        for (auto e : nodes) {
                SSNode n = dBG.getSSNode(e);
                int nodeMult = nodeCovModel.getExpMult(n.getAvgCov());

                int sumRight = 0;
                for (ArcIt r = n.rightBegin(); r != n.rightEnd(); r++) {
                        // count palindromic arcs double
                        int c = (r->getNodeID() == -e) ? 2 : 1;
                        sumRight += c * edgeCovModel.getExpMult(r->getCov());
                }

                if ((n.numRightArcs() > 0) && (nodeMult != sumRight))
                        return false;

                int sumLeft = 0;
                for (ArcIt l = n.leftBegin(); l != n.leftEnd(); l++) {
                        // count palindromic arcs double
                        int c = (l->getNodeID() == -e) ? 2 : 1;
                        sumLeft += c * edgeCovModel.getExpMult(l->getCov());
                }

                if ((n.numLeftArcs() > 0) && (nodeMult != sumLeft))
                        return false;
        }

        return true;
}

Multiplicity CRFSolver::getNodeMultiplicity(NodeRep node,
                                            const CovModel& nodeCovModel,
                                            const CovModel& edgeCovModel,
                                            int graphDepth)
{
        getSubgraph(node, graphDepth);

        // create a mapping between nodes/edges and CRF variables
        map<NodeRep, int> node2var;
        map<EdgeRep, int> edge2var;
        int varCRF = 0;
        for (const auto it : nodes){
                node2var[it] = varCRF++;
                //cout << varCRF -1 << "\t" << it.getNodeID() << endl;
        }
        for (const auto it : edges){
                edge2var[it] = varCRF++;
                //cout << varCRF -1 << "\t" << it.getSrcID() << "\t" << it.getDstID() << endl;
        }

        Multiplicity result;
        int targetCRFVar = node2var[node];
        if (solveSubgraph(node2var, edge2var, nodeCovModel,
                          edgeCovModel, result, targetCRFVar)) {
                perfCounter.totNumCRF++;
                perfCounter.totNumNodes += nodes.size();
                perfCounter.totNumEdges += edges.size();

                return result;
        }

        // fall-back to smaller subgraph if necessary
        return getNodeMultiplicity(node, nodeCovModel, edgeCovModel, graphDepth-1);
}

Multiplicity CRFSolver::getEdgeMultiplicity(EdgeRep edge,
                                            const CovModel& nodeCovModel,
                                            const CovModel& edgeCovModel,
                                            int graphDepth)
{
        getSubgraph(edge, graphDepth);
        
        // create a mapping between nodes/edges and CRF variables
        map<NodeRep, int> node2var;
        map<EdgeRep, int> edge2var;
        int varCRF = 0;
        for (const auto it : nodes)
                node2var[it] = varCRF++;
        for (const auto it : edges)
                edge2var[it] = varCRF++;
        
        Multiplicity result;
        int targetCRFVar = edge2var[edge];
        if (solveSubgraph(node2var, edge2var, nodeCovModel,
                edgeCovModel, result, targetCRFVar))
                return result;
        
        // fall-back to smaller subgraph if necessary
        return getEdgeMultiplicity(edge, nodeCovModel, edgeCovModel, graphDepth-1);
}

dai::FactorGraph CRFSolver::getLibdaiFG(const std::map<NodeRep, int>& node2var,
                                        const std::map<EdgeRep, int>& edge2var,
                                        std::vector<int>& firstMult,
                                        const CovModel& nodeCovModel,
                                        const CovModel& edgeCovModel)
{
        std::vector<dai::Factor> factorList;
        vector<int> card(node2var.size() + edge2var.size());
        int nodeSingletonF = 0;
        int arcSingletonF = 0;
        int flowF = 0;
        
        // first the nodes
        for (const auto it : node2var) {        
                SSNode node = dBG.getSSNode(it.first.getNodeID());
                const int& nodeVarID = it.second;
                
                int expMult = nodeCovModel.getExpMult(node.getAvgCov());
                firstMult[nodeVarID] = max(0, expMult - multMargin);
                card[nodeVarID] = expMult + multMargin - firstMult[nodeVarID] + 1;
                
                size_t numNbArcs = node.numLeftArcs() + node.numRightArcs();
                // add a singleton node factor in case:
                // - the node is bigger than 2k (less correlation with edges) OR
                // - there are no edges (in case of isolated nodes) OR
                // - there is a single incoming edge OR
                // - there is a single outgoing edge
                if ((node.getMarginalLength() < 2*Kmer::getK()) &&
                        (node.numLeftArcs() > 1) && (node.numRightArcs() > 1 && ! (edge2var.size() == 0)))
                        continue;
                
                int nIndptObs = max<int>(1, node.getMarginalLength() / (2*Kmer::getK()));
                
                dai::Factor F = createLDSingletonFactor(nodeVarID, card[nodeVarID],
                                                        firstMult[nodeVarID],
                                                        nodeCovModel, node.getAvgCov(),
                                                        nIndptObs);
                
                nodeSingletonF++;
                factorList.push_back(F);
        }
        
        for (const auto it : edge2var) {        // then the edges
                const EdgeRep& edge = it.first;
                const int& edgeVarID = it.second;
                
                double cov = dBG.getArc(edge).getCov();
                int expMult = edgeCovModel.getExpMult(cov);
                firstMult[edgeVarID] = max(0, expMult - multMargin);
                card[edgeVarID] = expMult + multMargin - firstMult[edgeVarID] + 1;
                
                // add a singleton edge factor in case:
                // - an adjacent node is not in the CRF (extremal edge) OR
                // - both adjacent nodes have multiple edges on that side
                if ((node2var.find(NodeRep(edge.getSrcID())) != node2var.end()) &&
                        (node2var.find(NodeRep(edge.getDstID())) != node2var.end()) && (! (node2var.size() == 0))) {
                        SSNode src = dBG.getSSNode(edge.getSrcID());
                if (src.numRightArcs() == 1)
                        continue;
                SSNode dst = dBG.getSSNode(edge.getDstID());
                if (dst.numLeftArcs() == 1)
                        continue;
                        }
                        
                        dai::Factor F = createLDSingletonFactor(edgeVarID, card[edgeVarID],
                                                                firstMult[edgeVarID],
                                                                edgeCovModel, cov, 1);
                        factorList.push_back(F);
                arcSingletonF++;
        }
        
        // create the flow conservation factors
        if ( !edge2var.empty() ){
                for (const auto it : node2var) {
                        NodeID currID = it.first.getNodeID();
                        const int& nodeVarID = it.second;
                        SSNode node = dBG.getSSNode(currID);
                        
                        // left flow conservation factor
                        vector<int> lEdgeVarID, lEdgeVarCard, lEdgeFirstMult;
                        vector<bool> lEdgePalindromic;
                        for (ArcIt lIt = node.leftBegin(); lIt != node.leftEnd(); lIt++) {
                                EdgeRep edge(lIt->getNodeID(), currID);
                                int edgeVarID = edge2var.at(edge);
                                
                                lEdgeVarID.push_back(edgeVarID);
                                lEdgeVarCard.push_back(card[edgeVarID]);
                                lEdgeFirstMult.push_back(firstMult[edgeVarID]);
                                lEdgePalindromic.push_back(lIt->getNodeID() == -currID);
                        }
                        
                        if (!lEdgeVarID.empty()) {
                                dai::Factor F = createLDFlowFactor(nodeVarID, card[nodeVarID],
                                                                   firstMult[nodeVarID], lEdgeVarID,
                                                                   lEdgeVarCard, lEdgeFirstMult,
                                                                   lEdgePalindromic, flowStrength);
                                factorList.push_back(F);
                                flowF++;
                        }
                        
                        
                        // right flow conservation factor
                        vector<int> rEdgeVarID, rEdgeVarCard, rEdgeFirstMult;
                        vector<bool> rEdgePalindromic;
                        for (ArcIt rIt = node.rightBegin(); rIt != node.rightEnd(); rIt++) {
                                EdgeRep edge(currID, rIt->getNodeID());
                                int edgeVarID = edge2var.at(edge);
                                
                                rEdgeVarID.push_back(edgeVarID);
                                rEdgeVarCard.push_back(card[edgeVarID]);
                                rEdgeFirstMult.push_back(firstMult[edgeVarID]);
                                rEdgePalindromic.push_back(currID == -rIt->getNodeID());
                        }
                        
                        if (!rEdgeVarID.empty()) {
                                dai::Factor F = createLDFlowFactor(nodeVarID, card[nodeVarID],
                                                                   firstMult[nodeVarID], rEdgeVarID,
                                                                   rEdgeVarCard, rEdgeFirstMult,
                                                                   rEdgePalindromic, flowStrength);
                                factorList.push_back(F);
                                flowF++;
                        }
                }
        }
        
        cout << "libDai Factor graph contains " << nodeSingletonF << " node singleton factors, "
        << arcSingletonF << " arc singleton factors, and " 
        << flowF << " flow factors" << endl;
        
        return dai::FactorGraph(factorList);
}

void CRFSolver::approxMult(NodeMap<Multiplicity>& nodeMult,
                        const CovModel& nodeCovModel,
                        EdgeMap<Multiplicity>& edgeMult,
                        const CovModel& edgeCovModel,
                        string modelID,
                        bool MAP, bool visualise)
{
        map<NodeRep, int> node2var;
        map<EdgeRep, int> edge2var;
        int varCRF = 0;
        for (const auto it : nodes){
                node2var[it] = varCRF++;
                //cout << varCRF -1 << "\t" << it.getNodeID() << endl;
        }
        for (const auto it : edges){
                edge2var[it] = varCRF++;
                //cout << varCRF -1 << "\t" << it.getSrcID() << "\t" << it.getDstID() << endl;
        }
        
        vector<int> firstMult(node2var.size() + edge2var.size());
        
        cout << "Getting libDAI FG representation ...." << endl;
        Util::startChrono();
        dai::FactorGraph fullFG = getLibdaiFG(node2var, edge2var, firstMult,
                                              nodeCovModel, edgeCovModel);
        cout << "... done (" << Util::stopChronoStr() << ")" << endl;
        
        
        // Properties for belief propagation TODO: lift this to settings or main.cpp and check appropriate values
        dai::PropertySet opts;
        if (Util::fileExists("libdai.props")) {
                ifstream options;
                options.open("libdai.props");
                string line;
                while(getline(options,line)) {
                        if (line.compare(0,1,"%") == 0)
                                continue;
                        size_t tab_pos = line.find("\t");
                        opts.set(line.substr(0,tab_pos),line.substr(tab_pos + 1, line.length()));
                }
        } else {
                opts.set("maxiter", (size_t) 500);
                opts.set("tol", 1e-3);
                opts.set("logdomain", true);
                opts.set("verbose", (size_t) 2);
                opts.set("maxtime", string("1800"));
                opts.set("updates", string("SEQMAX0L"));
                //opts.set("splashsize", (size_t) 5);
                opts.set("weightdecay", true);
                opts.set("forcestop", false);
                //opts.set("nthreads", (size_t) 1);
                //opts.set("resinit", string("UNIFORM"));
                //opts.set("damping", 0.2);
        }
        
        opts.set("modelcount", modelID);
        
        if (opts.getStringAs<size_t>("verbose") == 4 || visualise){
                map<size_t, size_t> nodevars;
                for (const auto it : nodes)
                        nodevars[node2var[it]] = it.getNodeID();
                cout << "writing FG representation ...." << endl;
                Util::startChrono();
                string fgfn = "factorGraph." + modelID + ".fg";
                fullFG.WriteToFile(fgfn.c_str(), 5);
                ofstream cytnodes, cytarcs;
                cytnodes.open("factorgraph." + modelID + ".cyt.nodes");
                cytarcs.open("factorgraph." + modelID + ".cyt.arcs");
                fullFG.dBGfg2Cytoscape(cytarcs, cytnodes,nodevars);
                cytnodes.close();
                cytarcs.close();
                cout << "... done (" << Util::stopChronoStr() << ")" << endl;
        }
        cout << "Starting libDAI LBP ...." << endl;
        Util::startChrono();
        double finalMaxDiff;
        vector<size_t> maxProbAssignment;
        map<int, dai::Factor> var2fact;
        map<int, int> label2idx;
        //if( opts.getStringAs<size_t>("nthreads") > 1 ){
        if ( ! opts.getStringAs<std::string>("updates").compare("SPLASH")){
                dai::ParallelBP bp;
                bp = dai::ParallelBP(fullFG, opts); 
                bp.init();
                finalMaxDiff = bp.run();
                
                // Get all beliefs by variable ID
                for (size_t i=0; i < fullFG.nrVars(); ++i) {
                        dai::Var v = fullFG.var(i);
                        int label = v.label();
                        dai::Factor belief = bp.beliefV(i);
                        var2fact[label] = belief;
                        label2idx[label] = i;
                }
        } else {
                dai::BP bp;
                if(MAP){
                        opts.set("updates", string("SEQMAX"));
                        bp = dai::BP(fullFG, opts("inference",string("MAXPROD")));
                } else {
                        bp = dai::BP(fullFG, opts); 
                }
                bp.init();
                finalMaxDiff = bp.run();
        /*if (finalMaxDiff > opts.getAs<double>("tol"))
                for (const auto it : nodes){
                        cout << varCRF -1 << "\t" << it.getNodeID() << endl;
                }
        */
                if(MAP)
                        maxProbAssignment = bp.findMaximum();
                
                // Get all beliefs by variable ID
                for (size_t i=0; i < fullFG.nrVars(); ++i) {
                        dai::Var v = fullFG.var(i);
                        int label = v.label();
                        dai::Factor belief = bp.beliefV(i);
                        var2fact[label] = belief;
                        label2idx[label] = i;
                }
        }
        cout << "... done (" << Util::stopChronoStr() << ")" << endl;
        
        // Extract multiplicities for the nodes
        for (int id = 0; id < nodes.size(); id ++) {
                NodeRep n = nodes[id];
                int label = node2var[n];
                dai::Factor f = var2fact[label];
                if(MAP){
                        Multiplicity mult(firstMult[label]  + maxProbAssignment[label2idx[label]]);
                        nodeMult[n] = mult;
                }else {
                        vector<double> logprob = f.log().p().p();
                        Multiplicity mult(firstMult[label], logprob);
                        nodeMult[n] = mult;
                }
                
        }
        
        // Extract multiplicities for the edges
        for (int id = 0; id < edges.size(); id ++) {
                EdgeRep n = edges[id];
                int label = edge2var[n];
                dai::Factor f = var2fact[label];
                if(MAP){
                        Multiplicity mult(firstMult[label]  + maxProbAssignment[label2idx[label]]);
                        edgeMult[n] = mult;
                }else {
                        vector<double> logprob = f.log().p().p();
                        Multiplicity mult(firstMult[label], logprob);
                        edgeMult[n] = mult;
                }
        }
}

void CRFSolver::approxSubgraphMult(NodeRep node,
                                   NodeMap<Multiplicity>& nodeMult,
                                   EdgeMap<Multiplicity>& edgeMult,
                                   const CovModel& nodeCovModel,
                                   const CovModel& edgeCovModel,
                                   int graphDepth, bool MAP)
{
	getSubgraph(node, graphDepth);

        string modelIdentifier = to_string(node.getNodeID()) + ".nb" + to_string(graphDepth);
	approxMult(nodeMult, nodeCovModel, edgeMult, edgeCovModel, modelIdentifier, MAP, true);
}

// ============================================================================
// ESTIMATE NODE/EDGE MULTIPLICITIES USING CONDITIONAL RANDOM FIELDS
// ============================================================================

void CRFMult::checkNodeFlow(CRFSolver& solver, WorkLoadBalancer& wlb,
                            const vector<NodeRep>& nodes,
                            const CovModel& nodeCovModel,
                            const CovModel& edgeCovModel,
                            vector<bool>& flowOK) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++)
                        flowOK[id] = solver.checkFlow(nodes[id],
                                                      nodeCovModel,
                                                      edgeCovModel,
                                                      maxGraphDepth);
}

void CRFMult::computeNodeMult(CRFSolver& solver,
                              WorkLoadBalancer& wlb,
                              const vector<NodeRep>& nodes,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              NodeMap<Multiplicity>& nodeMult) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++) {
                        NodeRep n = nodes[id];
                        nodeMult[n] = solver.getNodeMultiplicity(n,
                                nodeCovModel, edgeCovModel, maxGraphDepth);
                }
}

void CRFMult::computeEdgeMult(CRFSolver& solver,
                              WorkLoadBalancer& wlb,
                              const vector<EdgeRep>& edges,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              EdgeMap<Multiplicity>& edgeMult) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++) {
                        const EdgeRep& e = edges[id];
                        edgeMult[e] = solver.getEdgeMultiplicity(e,
                                nodeCovModel, edgeCovModel, maxGraphDepth);
                }
}

void CRFMult::checkFlow(const vector<NodeRep>& nodes,
                        vector<bool>& flowOK,
                        const CovModel& nodeCovModel,
                        const CovModel& edgeCovModel) const
{
        cout << "Checking flow for " << nodes.size() << " nodes (subgraph "
                "depth: " << maxGraphDepth << ")" << endl;
        vector<thread> wt(numThreads);
        vector<CRFSolver> solver(numThreads, CRFSolver(dBG, multMargin,
                                                       maxFactorSize, flowStrength));

        // assign a multiplicity to the nodes
        WorkLoadBalancer nodeWlb(0, nodes.size(), threadWork, "\tProcessing nodes");
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::checkNodeFlow, this, ref(solver[i]),
                               ref(nodeWlb), cref(nodes), cref(nodeCovModel),
                               cref(edgeCovModel), ref(flowOK));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));
}

void CRFMult::computeMult(NodeMap<Multiplicity>& nodeMult,
                          const CovModel& nodeCovModel,
                          EdgeMap<Multiplicity>& edgeMult,
                          const CovModel& edgeCovModel) const
{
        cout << "Computing multiplicity for " << nodeMult.size()
             << " nodes and " << edgeMult.size() << " edges (subgraph "
                "depth: " << maxGraphDepth << ")" << endl;
        vector<thread> wt(numThreads);
        vector<CRFSolver> solver(numThreads,
                CRFSolver(dBG, multMargin, maxFactorSize, flowStrength));

        // assign a multiplicity to the nodes
        vector<NodeRep> nodes; nodes.reserve(nodeMult.size());
        for (const auto& it : nodeMult)
                nodes.push_back(it.first);
        WorkLoadBalancer nodeWlb(0, nodes.size(), threadWork, "\tProcessing nodes");
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::computeNodeMult, this, ref(solver[i]),
                               ref(nodeWlb), cref(nodes), cref(nodeCovModel),
                               cref(edgeCovModel), ref(nodeMult));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));
        nodes.clear();

        // assign a multiplicity to the edges
        vector<EdgeRep> edges; edges.reserve(edgeMult.size());
        for (const auto& it : edgeMult)
                edges.push_back(it.first);
        WorkLoadBalancer edgeWlb(0, edges.size(), threadWork, "\tProcessing edges");
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::computeEdgeMult, this, ref(solver[i]),
                               ref(edgeWlb), cref(edges), cref(nodeCovModel),
                               cref(edgeCovModel), ref(edgeMult));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));

        // merge the performance counters
        for (const auto& it : solver)
                totPerfCounter.merge(it.perfCounter);
}

void CRFMult::approxMultAll(NodeMap<Multiplicity>& nodeMult,
                            const CovModel& nodeCovModel,
                            EdgeMap<Multiplicity>& edgeMult,
                            const CovModel& edgeCovModel,
                            bool MAP, bool singleCRF) const
{
        size_t ctr = 0;
        if (singleCRF) {
                vector<NodeRep> nodes = dBG.getNodeReps(dBG.getNumValidNodes());
                vector<EdgeRep> edges = dBG.getEdgeReps(dBG.getNumValidArcs());
                CRFSolver solver(dBG, multMargin, maxFactorSize, flowStrength, nodes, edges);
                solver.approxMult(nodeMult, nodeCovModel,
                                  edgeMult, edgeCovModel,
                                  "full." + to_string(ctr), MAP);
        } else {
                Bitvec handled(dBG.getNumNodes()+1);
                for(NodeRep nr : dBG.getNodeReps(dBG.getNumValidNodes())){
                        if (handled[abs(nr.getNodeID())])
                                continue;
                        CRFSolver solver(dBG, multMargin, maxFactorSize, flowStrength);
                        vector<NodeID> inSubgraph = solver.getFullConnectedComponent(nr);
                        if(inSubgraph.size() == 1){
                                nodeMult[nr] = solver.getNodeMultiplicity(nr, nodeCovModel, edgeCovModel, 0);
                                continue;
                        }
                        solver.approxMult(nodeMult, nodeCovModel,
                                        edgeMult, edgeCovModel,
                                        "full." + to_string(ctr), MAP);
                        for(NodeID id : inSubgraph)
                                handled[abs(id)] = true;
                        ctr++;
                }
        }
        
}

bool CRFMult::Estep(NodeMap<Multiplicity>& nodeMult,
                    const CovModel& nodeCovModel,
                    EdgeMap<Multiplicity>& edgeMult,
                    const CovModel& edgeCovModel,
                    double epsilon, bool approxInf,
                    bool map) const
{
        // store the old expected multiplicities to check for convergence
        NodeMap<int> oldNodeMult(nodeMult.size());
        for (const auto& it : nodeMult)
                oldNodeMult[it.first] = it.second.getExpMult();
        EdgeMap<int> oldEdgeMult(edgeMult.size());
        for (const auto& it : edgeMult)
                oldEdgeMult[it.first] = it.second.getExpMult();

        // compute (normalized) multiplicities given the model
        if (approxInf)
                approxMultAll(nodeMult, nodeCovModel,
                              edgeMult, edgeCovModel, map);
        else
                computeMult(nodeMult, nodeCovModel, edgeMult, edgeCovModel);

        // check for convergence
        size_t numNodeChange = 0;
        for (auto it : nodeMult)
                if (it.second.getExpMult() != oldNodeMult[it.first])
                        numNodeChange++;
        size_t numEdgeChange = 0;
        for (auto it : edgeMult)
                if (it.second.getExpMult() != oldEdgeMult[it.first])
                        numEdgeChange++;

                
        cout << "number of node changes: " << numNodeChange << endl;
        cout << "number of edge changes: " << numEdgeChange << endl;
        // convergence
        return (((double)numNodeChange/nodeMult.size() > epsilon) ||
                ((double)numEdgeChange/edgeMult.size() > epsilon));
}

bool CRFMult::Mstep(const NodeMap<Multiplicity>& nodeMult,
                    CovModel& nodeCovModel,
                    const EdgeMap<Multiplicity>& edgeMult,
                    CovModel& edgeCovModel, double epsilon, bool fixedZero) const
{
        bool retVal = true; // FIXME ?

        // ====================================================================
        // Compute the node model parameters
        // ====================================================================

        double nodeErrorLambda = nodeCovModel.getErrLambda();
        double nodeErrorODF = nodeCovModel.getErrorODF();
        double nodeErrorWeight = nodeCovModel.getWeight(0);
        
        if (! fixedZero){
        // We fit a negative binomial (NB) distribution to the error histogram.
        // The error histogram is always truncated: a) coverage zero is never
        // observed and b) low-coverage k-mers might have been removed by
        // a preprocessing tool like BCALM. We therefore use EM to fit the NB
        // parameters and infer the missing values from the spectrum.
        map<unsigned int, double> errorHist;
        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                if (it.second[0] > DOUBLE_SMALL) {
                        double f = node.getAvgCov() - floor(node.getAvgCov());
                        errorHist[node.getAvgCov()] += (1.0-f) * it.second[0];
                        errorHist[node.getAvgCov() + 1] += f * it.second[0];
                }
        }

        // We truncate the node histogram to the provided -abundance-min value.
        // This is relevant when using q-mer counts.
        unsigned int smallSupport = (dBG.getAbundanceMin() < 0)?
                errorHist.begin()->first : dBG.getAbundanceMin();
        for (unsigned int k = 0; k < smallSupport; k++)
                errorHist.erase(k);

        // The below EM procedure might not convergence, but the obtained NB
        // parameters should be sensible in all cases: only small support
        // values are inferred, therefore, the mean should never be larger than
        // the mean of the error histogram (and thus be close to zero).
        int maxIter = 10.0 / epsilon;
        int nIter = Util::fitTruncNegBinomEM(errorHist, nodeErrorLambda,
                                             nodeErrorODF, nodeErrorWeight,
                                             epsilon, maxIter);

        if (nIter > maxIter)
                cout << "\tWARNING: EM algorithm to fit node error model "
                        "did not converge\n";
        else
                cout << "\tEM algorithm to fit node error model converged "
                     << "after " << nIter << " iteration(s)" << endl;

        }
        // Compute node average and weights
        vector<double> wNode(nodeCovModel.getK(), 0.0);
        wNode[0] = nodeErrorWeight;

        double noml = 0.0, denl = 0.0;
        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                for (size_t m = 1; m < nodeCovModel.getK(); m++) {
                        wNode[m] += it.second[m];
                        noml += it.second[m] * node.getCov();
                        denl += it.second[m] * node.getMarginalLength() * m;
                }
        }

        double nodeLambda = max(0.5, noml / denl);

        // Compute node variance and ODF
        vector<double> wNodeSqr(nodeCovModel.getK(), 0.0);

        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                for (size_t m = 1; m < nodeCovModel.getK(); m++) {
                        double delta = node.getAvgCov() - m * nodeLambda;
                        wNodeSqr[m] += it.second[m] * delta * delta;
                }
        }
        double nom = 0.0, den = 0.0;
        for (size_t m = 1; m < nodeCovModel.getK(); m++) {
                nom += wNodeSqr[m] / (m * nodeLambda);
                den += wNode[m];
        }
        double nodeODF = nom / den;
        nodeODF = max(1.0, nodeODF);

        nodeCovModel = CovModel(nodeErrorLambda, nodeErrorODF,
                                nodeLambda, nodeODF, wNode);

        cout << "\tNode spectrum: " << nodeCovModel << endl;

        // ====================================================================
        // Compute the edge model parameters
        // ====================================================================

        double edgeErrorLambda = edgeCovModel.getErrLambda();
        double edgeErrorODF = edgeCovModel.getErrorODF();
        double edgeErrorWeight = edgeCovModel.getWeight(0);
        if (! fixedZero){
                map<unsigned int, double> errorHist;
                for (const auto& it : edgeMult) {
                        Arc& arc = dBG.getArc(it.first);
                        if (it.second[0] > DOUBLE_SMALL) {
                                double f = arc.getCov() - floor(arc.getCov());
                                errorHist[arc.getCov()] += (1.0-f) * it.second[0];
                                errorHist[arc.getCov() + 1] += f * it.second[0];
                        }
                }

                // We truncate the edge histogram on the same value as the node
                // histogram. Even when removing all k-mers with coverage < T, a few
                // arcs might still have a coverage < T. We do not want keep those
                // in our histogram as this will interfere with the EM algorithm.
                unsigned int smallSupport = (dBG.getAbundanceMin() < 0)?
                errorHist.begin()->first : dBG.getAbundanceMin();
                for (unsigned int k = 0; k < smallSupport; k++)
                        errorHist.erase(k);

                // The below EM procedure might not convergence, but the obtained NB
                // parameters should be sensible in all cases: only small support
                // values are inferred, therefore, the mean should never be larger than
                // the mean of the error histogram (and thus be close to zero).
                int maxIter = 10.0 / epsilon;
                int     nIter = Util::fitTruncNegBinomEM(errorHist, edgeErrorLambda,
                                                edgeErrorODF, edgeErrorWeight,
                                                epsilon, maxIter);

                if (nIter > maxIter)
                        cout << "\tWARNING: EM algorithm to fit edge error model "
                                "did not converge\n";
                else
                        cout << "\tEM algorithm to fit edge error model converged "
                        << "after " << nIter << " iteration(s)" << endl;
        }

        // Compute edge average and weights
        vector<double> wEdge(edgeCovModel.getK(), 0.0);
        wEdge[0] = edgeErrorWeight;

        vector<double> totWEdgeCov(edgeCovModel.getK(),0.0);
        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(it.first);
                for (size_t m = 1; m < edgeCovModel.getK(); m++) {
                        wEdge[m] += it.second[m];
                        totWEdgeCov[m] += it.second[m] * arc.getCov();
                }
        }

        noml = 0.0, denl = 0.0;
        for (size_t m = 1; m < nodeCovModel.getK(); m++) {
                noml += totWEdgeCov[m] / m;
                denl += wEdge[m];
        }

        double edgeLambda = max(0.5, noml / denl);

        // Compute edge variance and ODF
        vector<double> wEdgeSqr(edgeCovModel.getK(), 0.0);

        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(it.first);
                for (size_t m = 1; m < edgeCovModel.getK(); m++) {
                        double delta = arc.getCov() - m * edgeLambda;
                        wEdgeSqr[m] += it.second[m] * delta * delta;
                }
        }

        nom = 0.0, den = 0.0;
        for (size_t m = 1; m < edgeCovModel.getK(); m++) {
                nom += wEdgeSqr[m] / (m * edgeLambda);
                den += wEdge[m];
        }
        double edgeODF = nom / den;
        edgeODF = max(1.0, edgeODF);

        edgeCovModel = CovModel(edgeErrorLambda, edgeErrorODF,
                                edgeLambda, edgeODF, wEdge);

        cout << "\tEdge spectrum: " << edgeCovModel << endl;

        return retVal;
}

int CRFMult::computeMultEM(NodeMap<Multiplicity>& nodeMult,
                           CovModel& nodeCovModel,
                           EdgeMap<Multiplicity>& edgeMult,
                           CovModel& edgeCovModel,
                           double epsilon, int maxIter,
                           bool approxInf, bool map, bool fixedZero)
{
        int iter;
        for (iter = 1; iter <= maxIter; iter++) {
                cout << "Iteration " << iter << endl;
                // infer multiplicities given the model
                if (!Estep(nodeMult, nodeCovModel,
                           edgeMult, edgeCovModel,
                           epsilon, approxInf, map)) 
                        return iter;

                // update the model given the multiplicities
                if (!Mstep(nodeMult, nodeCovModel,
                           edgeMult, edgeCovModel, epsilon, fixedZero))
                                return iter;
        }

        return (iter);
}

