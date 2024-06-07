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

#include <sstream>
#include <random>
#include <utility>

#include "coverage.h"
#include "dbgraph.h"
#include "settings.h"
#include "util.h"

DSNode* SSNode::nodes = NULL;

using namespace std;

std::ostream &operator<<(std::ostream &out, const NodeChain& nc)
{
        for (size_t i = 0; i < nc.size(); i++)
                out << nc[i] << "\t";
        out << "(" << nc.getCount() << ")";
        return out;
}

bool DBGraph::consecutiveNPP(const NodePosPair& left,
                             const NodePosPair& right,
                             size_t offset) const
{
        // return false if one of the npps is invalid
        if (!left.isValid() || !right.isValid())
                return false;

        // left and right belong to the same node?
        if (left.getNodeID() == right.getNodeID())
                if ((left.getPosition() + offset) == right.getPosition())
                        return true;

        // left and right belong to connected nodes?
        SSNode lNode = getSSNode(left.getNodeID());
        if (lNode.rightArc(right.getNodeID()) == NULL)
                return false;

        // make sure the distance is consistent
        size_t dist = lNode.getMarginalLength() - left.getPosition() + right.getPosition();
        return (offset == dist);
}

void DBGraph::writeBCalm (const string& filename)
{
        ofstream ofs(filename.c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename
                     << " for writing" << endl;

        for (NodeID id = 1; id <= numNodes; id++)
        {
                SSNode node = getSSNode(id);
                if (! node.isValid())
                        continue;

                ofs << ">" << id << " LN:i:" <<
                node.getMarginalLength() + Kmer::getK() -1 << " KC:f:" <<
                node.getCov() << " km:f:" << node.getAvgCov();
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        ofs << "\tL:+:" << abs(rightID) << ":" << (rightID > 0 ? "+" : "-");
                }
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it ++) {
                        NodeID leftID = it->getNodeID();
                        ofs << "\tL:-:" << abs(leftID) << ":" << (leftID > 0 ? "-" : "+");
                }
                ofs << "\n";
                ofs << node.getSequence() << endl;
        }
        ofs.close();
}

void DBGraph::loadBCalm (const string& filename)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw ios_base::failure("Cannot open file " + filename);

        // figure out file size
        ifs.seekg(0, ios_base::end);
        size_t fileSize = ifs.tellg();
        ifs.clear();
        ifs.seekg(0, ios::beg);

        // first pass through the file to find out number of nodes and arcs
        Util::startChrono();
        string progressStr = "Reading file " + filename + " (first pass)";

        numNodes = numArcs = 0;

        string line;
        while (getline(ifs, line)) {
                if (numNodes % 1024 == 0)
                        Util::progress(progressStr, ifs.tellg(), fileSize);
                if (line.front() != '>')
                        continue;

                numNodes++;

                // every occurrence of "L:" denotes a new arc
                size_t pos = 0;
                while ((pos = line.find("L:", pos + 1, 2)) != string::npos)
                        numArcs++;
        }

        double elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);

        // +1 because index 0 is not used
        nodes = new DSNode[numNodes+1];
        SSNode::setNodePointer(nodes);

        // +2 because index 0 is not used, final index denotes 'end'.
        arcs = new Arc[numArcs+2];
        DSNode::setArcsPointer(arcs);

        // second pass through the file to store nodes and arcs
        cout << "Reading " << numNodes << " nodes and "
             << numArcs << " arcs..." << endl;

        ifs.clear();
        ifs.seekg(0, ios::beg);

        Util::startChrono();
        progressStr = "Reading file " + filename + " (second pass)";

        ArcID arcOffset = 1; NodeID nodeOffset = 1;
        while (getline(ifs, line)) {
                if (nodeOffset % 1024 == 0)
                        Util::progress(progressStr, ifs.tellg(), fileSize);
                DSNode& node = getDSNode(nodeOffset);

                // the line is a sequence
                if (line.front() != '>') {
                        node.setSequence(line);
                        nodeOffset++;
                        continue;
                }

                // the line is a FASTA descriptor line
                istringstream iss(line);

                // find "KC:i:" to figure out the k-mer coverage
                size_t pos = 0; int kmerCov = 0;
                if ((pos = line.find("KC:i:")) != string::npos) {
                        iss.seekg(pos + 5);
                        iss >> kmerCov;

                }

                node.setCov(kmerCov);

                vector<NodeID> leftArcs;
                vector<NodeID> rightArcs;

                // every occurrence of "L:" denotes a new arc
                pos = 0;
                while ((pos = line.find("L:", pos + 1, 2)) != string::npos) {
                        iss.seekg(pos + 2);
                        char c, l, r;
                        int dstID;
                        iss >> l >> c >> dstID >> c >> r >> right;
                        dstID++;        // we number from 1, not 0

                        if (l == '+')
                                rightArcs.push_back(r == '+' ? dstID : -dstID);
                        else
                                leftArcs.push_back(r == '-' ? dstID : -dstID);
                }

                node.setNumLeftArcs(leftArcs.size());
                node.setFirstLeftArcID(arcOffset);
                for (auto dstID : leftArcs)
                        arcs[arcOffset++].setNodeID(dstID);

                node.setNumRightArcs(rightArcs.size());
                node.setFirstRightArcID(arcOffset);
                for (auto dstID : rightArcs)
                        arcs[arcOffset++].setNodeID(dstID);
        }

        elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);
        ifs.close();

        // figure out the value for k
        bool autoDetectK = false;

        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;

                SSNode node = getSSNode(id);

                if (node.numRightArcs() < 2)
                        continue;

                ArcIt it = node.rightBegin();
                SSNode nb1 = getSSNode(it->getNodeID());
                it++;
                SSNode nb2 = getSSNode(it->getNodeID());
                string s1 = nb1.getSequence();
                string s2 = nb2.getSequence();

                int k = 0;
                while (s1[k] == s2[k])          // overlap == k-1
                        k++;
                k++;                            // add 1 to get k

                Kmer::setWordSize(k);
                autoDetectK = true;
                break;
        }

        if (!autoDetectK)
                throw runtime_error("Cannot infer kmer size from input file");
        else
                cout << "Kmer size is " << Kmer::getK() << endl;

        numValidNodes = numNodes;
        numValidArcs = numArcs;
}

void DBGraph::writeBinary(const std::string& filename) const
{
        ofstream ofs(filename.c_str(), ios::binary);
        // write k
        size_t k = Kmer::getK();
        ofs.write((char*)&k, sizeof(k));

        // A) write nodes
        ofs.write((char*)&numNodes, sizeof(numNodes));
        for (NodeID id = 1; id <= numNodes; id++)
                getDSNode(id).write(ofs);

        // B) write arcs
        ofs.write((char*)&numArcs, sizeof(numArcs));
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].write(ofs);

        ofs.close();

        cout << "Wrote " << numNodes << " nodes and "
             << numArcs << " arcs" << endl;
}

void DBGraph::loadBinary(const std::string& filename)
{
        // read the metadata
        ifstream ifs(filename.c_str(), ios::binary);
        if (!ifs)
                throw ios_base::failure("Cannot open " + filename);

        size_t k;
        ifs.read((char*)&k, sizeof(k));
        Kmer::setWordSize(k);

        // A) create the nodes
        ifs.read((char*)&numNodes, sizeof(numNodes));
        // +1 because index 0 isn't used
        nodes = new DSNode[numNodes+1];
        SSNode::setNodePointer(nodes);
        for (NodeID id = 1; id <= numNodes; id++)
                getDSNode(id).read(ifs);

        // B) create the arcs
        ifs.read((char*)&numArcs, sizeof(numArcs));
        // +2 because index 0 isn't used, final index denotes 'end'
        arcs = new Arc[numArcs+2];
        DSNode::setArcsPointer(arcs);
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].read(ifs);

        ifs.close();

        numValidNodes = numNodes;
        numValidArcs = numArcs;
}

void DBGraph::writeContigs(const std::string& filename) const
{
        ofstream ofs(filename.c_str());

        size_t contigID = 0;
        for (NodeID id = 1; id <= getNumNodes(); id++) {
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;

                //ofs << ">contig_" << contigID++ << "\n";
                ofs << ">contig_" << id << "\n";
                Util::writeSeqWrap(ofs, n.getSequence(), 60);
        }
}

void DBGraph::getGraph(std::vector<NodeID>& nodes, std::vector<EdgeID>& edges)
{
        nodes.reserve(2 * getNumValidNodes());
        edges.reserve(getNumValidArcs());

        for (NodeID srcID = -numNodes; srcID <= numNodes; srcID++) {
                if (srcID == 0)
                        continue;
                SSNode n = getSSNode(srcID);
                if (!getSSNode(srcID).isValid())
                        continue;
                nodes.push_back(srcID);

                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID dstID = it->getNodeID();
                        edges.push_back(make_pair(srcID, dstID));
                }
        }
}

void DBGraph::getSubgraph(NodeID seedNode, vector<NodeID>& nodes,
                          vector<EdgeID>& edges, size_t maxDepth) const
{
        priority_queue<NodeDFS, vector<NodeDFS>, NodeDFSComp> todo;
        todo.push(NodeDFS(seedNode, 0));
        set<NodeID> nodeSet;

        while (!todo.empty()) {
                // get and erase the current node
                NodeDFS currTop = todo.top();
                todo.pop();
                NodeID thisID = currTop.nodeID;
                size_t thisDepth = currTop.depth;

                SSNode n = getSSNode(thisID);

                // if the node was already handled, skip
                if (nodeSet.find(thisID) != nodeSet.end())
                        continue;

                nodeSet.insert(thisID);

                // don't go any deeper if we've reached the maximum depth
                if (thisDepth >= maxDepth)
                        continue;

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        if (nodeSet.find(rightID) != nodeSet.end())
                                continue;       // edge already added by rightID

                        edges.push_back(make_pair(thisID, rightID));
                        todo.push(NodeDFS(rightID, thisDepth+1));
                }

                // process the left arcs
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        if (nodeSet.find(leftID) != nodeSet.end())
                                continue;       // edge already added by leftID

                        edges.push_back(make_pair(leftID, thisID));
                        todo.push(NodeDFS(leftID, thisDepth + 1));
                }
                
        }
        
        nodes = vector<NodeID>(nodeSet.begin(), nodeSet.end());
}

/*void DBGraph::getSubgraph(priority_queue<NodeRepDepth, vector<NodeRepDepth>,
                          NodeRepComp>& todo, set<NodeRep>& nodes,
                          set<EdgeRep>& edges, size_t maxDepth) const
{
        while (!todo.empty()) {
                // get and erase the current node
                NodeRepDepth currTop = todo.top();
                todo.pop();
                NodeRep thisID = currTop.nodeRep;
                size_t thisDepth = currTop.depth;

                // if the node was already handled, skip
                if (nodes.find(thisID) != nodes.end())
                        continue;

                // mark this node as handled
                nodes.insert(thisID);

                // process the right arcs
                SSNode node = getSSNode(thisID);
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        edges.insert(EdgeRep(thisID, rightID));
                        if (thisDepth < maxDepth)
                                todo.push(NodeRepDepth(rightID, thisDepth+1));
                }

                // process the left arcs
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        edges.insert(EdgeRep(leftID, thisID));
                        if (thisDepth < maxDepth)
                                todo.push(NodeRepDepth(leftID, thisDepth + 1));
                }
        }
}

void DBGraph::getSubgraph(NodeRep seedNode, set<NodeRep>& nodes,
                          set<EdgeRep>& edges, size_t maxDepth) const
{
        // a queue containing nodeIDs to handle + their depth
        priority_queue<NodeRepDepth, vector<NodeRepDepth>, NodeRepComp> todo;
        todo.push(NodeRepDepth(seedNode, 0));
        
        getSubgraph(todo, nodes, edges, maxDepth);
}
                          
void DBGraph::getSubgraph(EdgeRep seedEdge, set<NodeRep>& nodes,
                        set<EdgeRep>& edges, size_t maxDepth) const
{
        edges.insert(seedEdge);
        if (maxDepth == 0)
                return;
        
        // a queue containing nodeIDs to handle + their depth
        priority_queue<NodeRepDepth, vector<NodeRepDepth>, NodeRepComp> todo;
        todo.push(NodeRepDepth(seedEdge.getSrcID(), 1));
        todo.push(NodeRepDepth(seedEdge.getDstID(), 1));
        
        getSubgraph(todo, nodes, edges, maxDepth);
}*/


vector<NodeRep> DBGraph::getNodeReps(const string& filename) const
{
        ifstream ifs(filename.c_str());

        NodeID n;
        size_t numNodes = 0;
        while (ifs >> n)
                numNodes++;

        vector<NodeRep> nodes;
        nodes.reserve(numNodes);

        ifs.clear();
        ifs.seekg(0, ios::beg);

        while (ifs >> n) {
                if (!nodeExists(n))
                        throw runtime_error("Node with ID " + to_string(n) + " does not exist");
                nodes.push_back(NodeRep(n));
        }

        return nodes;
}

vector<NodeRep> DBGraph::getNodeReps(size_t N) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++)
                if (getSSNode(i).isValid())
                        nodes.push_back(NodeRep(i));

        N = min(N, nodes.size());
        if (N == nodes.size())  // if you need all nodes we do not shuffle
                return nodes;

        // sample N nodes using the Fisher-Yates algoritm
        random_device rd;
        mt19937 mt(rd());
        for (size_t i = 0; i < N; i++) {
                uniform_int_distribution<size_t> dis(i, nodes.size() - 1);
                swap(nodes[i], nodes[dis(mt)]);
        }

        return vector<NodeRep>(nodes.begin(), nodes.begin() + N);
}

vector<EdgeRep> DBGraph::getEdgeReps(const string& filename) const
{
        ifstream ifs(filename.c_str());

        ArcID l, r;
        size_t numEdges = 0;
        while (ifs >> l >> r)
                numEdges++;

        vector<EdgeRep> edges;
        edges.reserve(numEdges);

        ifs.clear();
        ifs.seekg(0, ios::beg);

        while (ifs >> l >> r) {
                if (!edgeExists(EdgeRep(l, r)))
                        throw runtime_error("Edge from ID " + to_string(l) + " to ID " + to_string(r) + " does not exist");
                edges.push_back(EdgeRep(l, r));
        }

        return edges;
}

vector<EdgeRep> DBGraph::getEdgeReps(size_t N) const
{
        vector<EdgeRep> edges;
        edges.reserve(getNumValidArcs() / 2);   // roughly half of the nodes
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        if (id <= -it->getNodeID())
                                edges.push_back(EdgeRep(id, it->getNodeID()));
        }

        N = min(N, edges.size());
        if (N == edges.size())  // if you need all edges we do not shuffle
                return edges;

        // sample N edges using the Fisher-Yates algoritm
        random_device rd;
        mt19937 mt(rd());
        for (size_t i = 0; i < N; i++) {
                uniform_int_distribution<size_t> dis(i, edges.size() - 1);
                swap(edges[i], edges[dis(mt)]);
        }

        return vector<EdgeRep>(edges.begin(), edges.begin() + N);
}

vector<EdgeRep> DBGraph::getLowCovEdges(double threshold) const
{
        vector<EdgeRep> edges;
        edges.reserve(getNumValidArcs() / 2);   // roughly half of the nodes
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        if (id <= -it->getNodeID() && it->getCov() <= threshold)
                                edges.push_back(EdgeRep(id, it->getNodeID()));
                }
        }

        return edges;
}

vector<NodeRep> DBGraph::getLowCovNodes(double threshold) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;

                if (n.getAvgCov() <= threshold)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

vector<NodeRep> DBGraph::getLowCovTips(double threshold, size_t maxLen) const
{
        vector<NodeRep> nodes;
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;       // deleted node

                if ((n.numLeftArcs() > 0) && (n.numRightArcs() > 0))
                        continue;       // not a tip

                if ((maxLen != 0) && (n.getMarginalLength() > maxLen))
                        continue;       // node too long

                if (n.getAvgCov() <= threshold)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

vector<NodeRep> DBGraph::getLowCovBubbles(double threshold) const
{
        vector<NodeRep> nodes;
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;       // deleted node

                if ((n.numLeftArcs() != 1) || (n.numRightArcs() != 1))
                        continue;       // not a bubble

                /*NodeID leftID = n.leftBegin()->getNodeID();
                NodeID rightID = n.rightBegin()->getNodeID();

                bool isBubble = false;  // find parallel path
                SSNode l = getSSNode(leftID);
                for (ArcIt it = l.rightBegin(); it != l.rightEnd(); it++) {
                        if (it->getNodeID() == i)
                                continue;
                        if (edgeExists(EdgeRep(it->getNodeID(), rightID)))
                                isBubble = true;
                }

                if (!isBubble)
                        continue;*/

                if (n.getMarginalLength() != Kmer::getK())
                        continue;       // not a bubble

                if (n.getAvgCov() <= threshold)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

void DBGraph::getNodeCovHist(map<int, double>& hist)
{
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;

                hist[node.getAvgCov() + 0.5]++;// node.getMarginalLength();
        }
}

void DBGraph::writeCytoscapeGraph(const std::string& filename,
                                  vector<NodeID> nodes,
                                  vector<pair<NodeID, NodeID> > edges,
                                  const NodeMap<Multiplicity>& enm,
                                  const EdgeMap<Multiplicity>& eem,
                                  const NodeMap<int>& tnm,
                                  const EdgeMap<int>& tem) const
{
        // A) write all arcs
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs"
                     << " for writing" << endl;
        ofs << "Source node\tTarget node\tCoverage\tTrue mult.\t"
               "Est mult.\tLOR\tBelief\tcorrect\n";

        for (const auto& edge : edges) {
                EdgeRep er(edge);
                ArcID id = getArcID(er);

                auto teIt = tem.find(er);
                int trueMult = (teIt == tem.end()) ? -1 : teIt->second;

                auto eeIt = eem.find(er);
                int estMult = (eeIt == eem.end()) ? -1 : eeIt->second.getExpMult();

                double multLOR = (eeIt == eem.end()) ?
                        -1.0 : eeIt->second.getExpMultLogOR();
                double multBelief = (eeIt == eem.end()) ?
                        -1.0 : exp(eeIt->second.getExpMultLProb());

                ofs << edge.first << "\t" << edge.second << "\t"
                    << getArc(id).getCov() << "\t" << trueMult << "\t"
                    << estMult << "\t" << multLOR << "\t" << multBelief << "\t" << (trueMult==estMult? 1:0)<< "\n";
        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes"
                     << " for writing" << endl;

        ofs << "Node ID\tMarginal length\tLeft arcs\tRight arcs\tCoverage\t"
               "True mult.\tEst. mult.\tLOR\tBelief\tcorrect\tSequence\n";

        for (const auto& id : nodes) {
                NodeRep nr(id);
                SSNode n = getSSNode(id);

                auto tnIt = tnm.find(nr);
                int trueMult = (tnIt == tnm.end()) ? -1 : tnIt->second;

                auto enIt = enm.find(nr);
                int estMult = (enIt == enm.end()) ? -1 : enIt->second.getExpMult();

                double multLOR = (enIt == enm.end()) ?
                        -1.0 : enIt->second.getExpMultLogOR();
                double multBelief = (enIt == enm.end()) ?
                        -1.0 : exp(enIt->second.getExpMultLProb());

                ofs << id << "\t" << n.getMarginalLength() << "\t"
                    << (int)n.numLeftArcs() << "\t" << (int)n.numRightArcs()
                    << "\t" << n.getAvgCov() << "\t" << trueMult << "\t"
                    << estMult << "\t" << multLOR << "\t" << multBelief<< "\t" << (trueMult==estMult? 1:0) << "\t"
                    << n.getSequence() << "\n";
        }
        ofs.close();
}

void DBGraph::defragNodes()
{
        vector<NodeID> old2new(numNodes + 1, 0);

        // defragment the node array
        for (NodeID oldID = 1, newID = 1; oldID <= numNodes; oldID++) {
                SSNode n = getSSNode(oldID);
                if (!n.isValid())
                        continue;

                old2new[oldID] = newID;

                // we use the move assignment operator for efficiency
                // (to avoid a deep copy of the TString member)
                nodes[newID] = move(nodes[oldID]);
                newID++;
        }

        // update the arcs to point to the new nodeIDs
        for (ArcID id = 1; id <= numArcs; id++) {
                NodeID oldID = arcs[id].getNodeID();
                NodeID newID = (oldID < 0) ? -old2new[-oldID] : old2new[oldID];
                arcs[id].setNodeID(newID);
        }

        numNodes = numValidNodes;
        cout << "Defragged nodes array: " << numNodes << " nodes\n";
}

void DBGraph::defragArcs()
{
        vector<ArcID> old2new(numArcs + 2, 0);

        numValidArcs = 0;
        for (ArcID oldID = 1; oldID <= numArcs + 1; oldID++) {
                if (!arcs[oldID].isValid())
                        continue;

                numValidArcs++;
                old2new[oldID] = numValidArcs;

                // we use the regular assignment operator
                arcs[numValidArcs] = arcs[oldID];
        }

        // update the nodes to point to the new arcIDs
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;

                ArcID oldID = n.getFirstLeftArcID();
                ArcID newID = (oldID < 0) ? -old2new[-oldID] : old2new[oldID];
                n.setFirstLeftArcID(newID);

                oldID = n.getFirstRightArcID();
                newID = (oldID < 0) ? -old2new[-oldID] : old2new[oldID];
                n.setFirstRightArcID(newID);
        }

        numArcs = numValidArcs;
        cout << "Defragged arcs array: " << numArcs << " arcs\n";
}

void DBGraph::createKmerNPPTable(KmerNPPTable& table) const
{
        Util::startChrono();
        string progressStr("Populating <k-mer, node> table");

        // count the number of k-mers in the graph
        size_t numKmers = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode& node = nodes[id];
                if (!node.isValid())
                        continue;
                numKmers += node.getMarginalLength();
        }

        table.resize(numKmers);

        // populate the table with kmers
        for (NodeID id = 1; id <= numNodes; id++) {
                if (id % 1024 == 0)
                        Util::progress(progressStr, id, numNodes);

                const DSNode &node = nodes[id];
                if (!node.isValid())
                        continue;
                const TString& tStr = node.getTSequence();
                Kmer kmer(tStr);
                NodePosition pos = 0;
                table.insert(kmer, NodePosPair(id, pos++));

                for (size_t i = Kmer::getK(); i < tStr.getLength(); i++) {
                        kmer.pushNucleotideRight(tStr[i]);
                        table.insert(kmer, NodePosPair(id, pos++));
                }
        }

        double elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);
}

int DBGraph::getAbundanceMin() const {
        return settings.getAbundanceMin();
}
