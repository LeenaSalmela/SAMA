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

#ifndef SSNODE_H
#define SSNODE_H

#include "kmer/tkmer.h"
#include "dsnode.h"

#include <vector>
#include <google/sparse_hash_map>

// ============================================================================
// NODE REPRESENTATIVE CLASS
// ============================================================================

/**
 * Represent a node in a graph as a nodeID. Takes into account that identifiers
 * nodeID/-nodeID represent the same physical nodes
 */
class NodeRep {

private:
        NodeID nodeID;          // node identifier

public:
        /**
         * Default constructor
         * @param nodeID_ node identifier
         */
        NodeRep(NodeID nodeID_) : nodeID(abs(nodeID_)) {}

        NodeRep(){}

        /**
         * Get the node identifier
         * @return the node identifier
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Operator == overloading
         * @return true or false
         */
        bool operator==(const NodeRep& rhs) const {
                return (nodeID == rhs.nodeID);
        }

        /**
         * Operator < overloading
         * @return true or false
         */
        bool operator<(const NodeRep& rhs) const {
                return nodeID < rhs.nodeID;
        }

        /**
         * Implicit conversion back to NodeID
         * @return Node identifier
         */
        operator NodeID() const {
                return nodeID;
        }

        /**
         * Compute a hash function
         * @return Hash value
         */
        size_t getHash() const {

                size_t w = nodeID;
                size_t hash = 0;
                w = ~w + (w << 21);             // key = (key << 21) - key - 1;
                w = w ^ (w >> 24);
                w = (w + (w << 3)) + (w << 8);  // key * 265
                w = w ^ (w >> 14);
                w = (w + (w << 2)) + (w << 4);  // key * 21
                w = w ^ (w >> 28);
                w = w + (w << 31);
                hash = hash ^ size_t(w);

                return hash;
        }
};

// ============================================================================
// HASH FUNCTION
// ============================================================================

struct NodeHash {
        size_t operator()(const NodeRep &nr) const {
                return nr.getHash();
        }
};

// ============================================================================
// NODE MAP
// ============================================================================

template<class T>
using NodeMap = google::sparse_hash_map<NodeRep, T, NodeHash>;

// ============================================================================
// SINGLE STRANDED NODE CLASS
// ============================================================================

class SSNode {

private:
        static DSNode *nodes;   // pointer to the double stranded nodes

        NodeID nodeID;          // identiffier of the node
        DSNode *dsNode;         // reference to the double stranded node

public:
        /**
         * Default constructor
         */
        SSNode() : nodeID(0), dsNode(NULL) {}

        /**
         * Constructor
         * @param dsNode Double stranded node
         * @param ID Unique identifier of the node
         */
        SSNode(DSNode* dsNode, NodeID nodeID) : nodeID(nodeID), dsNode(dsNode) {
                assert(nodeID != 0);
        }

        /**
         * Constructor
         * @param it Iterator pointing to this node
         */
        SSNode(const ArcIt& it) : nodeID(it->getNodeID()), dsNode(nodes + abs(nodeID)) {
                assert(nodeID != 0);
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag1(bool flag) {
                dsNode->setFlag1(flag);
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag1() const {
                return dsNode->getFlag1();
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag2(bool flag) {
                dsNode->setFlag2(flag);
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag2() const {
                return dsNode->getFlag2();
        }

#ifndef AVG_COV
        /**
         * Set the count
	 * @param pos The position
         * @param target The coverage
         */
        void setCount(int pos, Count target) {
	  if (nodeID > 0) {
	    dsNode->setCount(pos, target);
	  } else {
	    dsNode->setCount(getMarginalLength()-pos-1, target);
	  }
        }

        void clearCounts() {
	  dsNode->clearCounts();
	}

  
        /**
         * Get the count
	 * @param pos The position
         * @return The coverage
         */
        Coverage getCount(int pos) const {
	  if (nodeID > 0) {
	    return dsNode->getCount(pos);
	  } else {
	    return dsNode->getCount(getMarginalLength()-pos-1);
	  }
        }

        /**
         * Atomically increment the count
	 * @param pos The position
         * @param rhs Right hand side (default = 1.0f)
         */
        void incCount(int pos, Coverage rhs = 1.0f) {
	  if (nodeID > 0) {
	    dsNode->incCount(pos, rhs);
	  } else {
	    dsNode->incCount(getMarginalLength()-pos-1, rhs);
	  }
        }

        /**
         * Set the arc count
	 * @param pos The position
         * @param target The coverage
         */
        void setArcCount(int pos, Count target) {
	  if (nodeID > 0) {
	    dsNode->setArcCount(pos, target);
	  } else {
	    dsNode->setArcCount(getMarginalLength()-pos-2, target);
	  }
        }

        void clearArcCounts() {
	  dsNode->clearArcCounts();
	}

  
        /**
         * Get the arc count
	 * @param pos The position
         * @return The coverage
         */
        Coverage getArcCount(int pos) const {
	  if (nodeID > 0) {
	    return dsNode->getArcCount(pos);
	  } else {
	    return dsNode->getArcCount(getMarginalLength()-pos-2);
	  }
        }

        /**
         * Atomically increment the arc count
	 * @param pos The position
         * @param rhs Right hand side (default = 1.0f)
         */
        void incArcCount(int pos, Coverage rhs = 1.0f) {
	  if (nodeID > 0) {
	    dsNode->incArcCount(pos, rhs);
	  } else {
	    dsNode->incArcCount(getMarginalLength()-pos-2, rhs);
	  }
        }

        void initializeCounts() {
	  dsNode->initializeCounts();
	}
  
#endif

        /**
         * Set the coverage
         * @param target The coverage
         */
        void setCov(Coverage target) {
                dsNode->setCov(target);
        }

        /**
         * Get the coverage
         * @return The coverage
         */
        Coverage getCov() const {
                return dsNode->getCov();
        }

        /**
         * Atomically increment the coverage
         * @param rhs Right hand side (default = 1.0f)
         */
        void incCov(Coverage rhs = 1.0f) {
                dsNode->incCov(rhs);
        }

        /**
         * Get the avarge kmer coverage
         * @return The average kmer coverage
         */
        double getAvgCov() const {
                return (double)dsNode->getCov() / (double)dsNode->getMarginalLength();
        }

        /**
         * Invalidate this node
         */
        void invalidate() {
                dsNode->invalidate();
        }

        /**
         * Check if a node is invalidated
         * @return True or false
         */
        bool isValid() const {
                if (nodeID == 0)
                        return false;
                return dsNode->isValid();
        }

        /**
         * Get the identifier of this node
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Operator '==' overloading
         * @param rhs Right hand side SSNode
         * @return True if they're equal
         */
        bool operator==(const SSNode &rhs) const {
                if (dsNode != rhs.dsNode)
                        return false;
                return (nodeID == rhs.nodeID);
        }

        /**
         * Operator '!=' overloading
         * @param rhs Right hand side kmer
         * @return True if they're different
         */
        bool operator!=(const SSNode &rhs) const {
                return !(*this == rhs);
        }

        /**
         * Get the length of the node (# nucleotides in DNA string)
         * @return The length of the node
         */
        size_t length() const {
                return dsNode->getLength();
        }

        /**
         * The marginal length == length - k + 1
         * @return The marginal length of the node
         */
        size_t getMarginalLength() const {
                return dsNode->getMarginalLength();
        }

        /**
         * Get the number of left arcs
         * @return The number of left arcs
         */
        uint8_t numLeftArcs() const {
                return (nodeID > 0) ?
                        dsNode->numLeftArcs() : dsNode->numRightArcs();
        }

        /**
         * Get the number of right arcs
         * @return The number of right arcs
         */
        uint8_t numRightArcs() const {
                return (nodeID > 0) ?
                        dsNode->numRightArcs() : dsNode->numLeftArcs();
        }

        /**
         * Get the number of left arcs
         * @param The number of left arcs
         */
        void setNumLeftArcs(uint8_t numArcs) const {
                if (nodeID > 0)
                        dsNode->setNumLeftArcs(numArcs);
                else
                        dsNode->setNumRightArcs(numArcs);
        }

        /**
         * Get the number of right arcs
         * @param The number of right arcs
         */
        void setNumRightArcs(uint8_t numArcs) const {
                if (nodeID > 0)
                        dsNode->setNumRightArcs(numArcs);
                else
                        dsNode->setNumLeftArcs(numArcs);
        }

        void swapRightArcsSign() {
                if (nodeID > 0)
                        dsNode->swapRightArcsSign();
                else
                        dsNode->swapLeftArcsSign();
        }

        void copyRightArcs(SSNode &source) {
                int numRightArcs = source.numRightArcs();
                source.setNumRightArcs(0);
                setNumRightArcs(numRightArcs);
                setFirstRightArcID(source.getFirstRightArcID());
                // if they have opposite sign
                if ((nodeID < 0) != (source.getNodeID() < 0))
                        swapRightArcsSign();
        }

        /**
         * Get an iterator pointing to the first left arc
         * @return an iterator pointing to the first left arc
         */
        ArcIt leftBegin() const {
                return (nodeID > 0) ?
                        dsNode->leftBegin(false) : dsNode->rightBegin(true);
        }

        /**
         * Get an iterator pointing past the last left arc
         * @retur nan iterator pointing to the last left arc
         */
        ArcIt leftEnd() const {
                return (nodeID > 0) ?
                        dsNode->leftEnd(false) : dsNode->rightEnd(true);
        }

        /**
         * Get an iterator pointing to the first right arc
         * @return an iterator pointing to the first left arc
         */
        ArcIt rightBegin() const {
                return (nodeID > 0) ?
                        dsNode->rightBegin(false) : dsNode->leftBegin(true);
        }

        /**
         * Get an iterator pointing past the last right arc
         * @retur nan iterator pointing to the last right arc
         */
        ArcIt rightEnd() const {
                return (nodeID > 0) ?
                        dsNode->rightEnd(false) : dsNode->leftEnd(true);
        }

        /**
         * Delete the left arcs
         */
        void deleteAllLeftArcs() {
                if (nodeID > 0)
                        dsNode->deleteLeftArcs();
                else
                        dsNode->deleteRightArcs();
        }

        /**
         * Delete the right arcs
         */
        void deleteAllRightArcs() {
                if (nodeID > 0)
                        dsNode->deleteRightArcs();
                else
                        dsNode->deleteLeftArcs();
        }

        /**
         * Get the identifier for the left right arc
         * @return The identifier for the left right arc
         */
        ArcID getFirstLeftArcID() {
                if (nodeID > 0)
                        return dsNode->getFirstLeftArcID();
                return dsNode->getFirstRightArcID();
        }

        /**
         * Get the identifier for the first right arc
         * @return The identifier for the first right arc
         */
        ArcID getFirstRightArcID() {
                if (nodeID > 0)
                        return dsNode->getFirstRightArcID();
                return dsNode->getFirstLeftArcID();
        }

        /**
         * Set the identifier for the left right arc
         * @param The identifier for the left right arc
         */
        void setFirstLeftArcID(ArcID target) {
                if (nodeID > 0)
                        return dsNode->setFirstLeftArcID(target);
                return dsNode->setFirstRightArcID(target);
        }

        /**
         * Set the identifier for the first right arc
         * @param The identifier for the first right arc
         */
        void setFirstRightArcID(ArcID target) {
                if (nodeID > 0)
                        return dsNode->setFirstRightArcID(target);
                return dsNode->setFirstLeftArcID(target);
        }

        /**
         * Delete a specific left arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteLeftArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->deleteLeftArc(targetID);
                return dsNode->deleteRightArc(-targetID);
        }

        /**
         * Delete a specific right arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteRightArc(NodeID targetID) {
                if (nodeID > 0)
                        return dsNode->deleteRightArc(targetID);
                return dsNode->deleteLeftArc(-targetID);
        }

        /**
         * Get a specific left arc
         * @param targetID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* leftArc(NodeID targetID) const {
                if (nodeID > 0)
                        return dsNode->leftArc(targetID);
                return dsNode->rightArc(-targetID);
        }

        /**
         * Get a specific right arc
         * @param targetID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* rightArc(NodeID targetID) const {
                if (nodeID > 0)
                        return dsNode->rightArc(targetID);
                return dsNode->leftArc(-targetID);
        }

        /**
         * Get the nodeID a specific left arc
         * @param nucleotide Nucleotide under consideration
         * @return Node identifier
         */
       NodeID getLeftArcNodeID(char nucleotide) const {
                for (ArcIt it = leftBegin(); it != leftEnd(); it++) {
                        if (SSNode(it).peekNucleotideMarginalRight() == nucleotide)
                                return it->getNodeID();
                }

                return 0;
        }

        /**
         * Get a specific right arc
         * @param nucleotide Nucleotide under consideration
         * @return Pointer to the specific arc, NULL if not found
         */
        NodeID getRightArcNodeID(char nucleotide) const {
                for (ArcIt it = rightBegin(); it != rightEnd(); it++)
                        if (SSNode(it).peekNucleotideMarginalLeft() == nucleotide)
                                return it->getNodeID();

                return 0;
        }

        /**
         * Replace a left arc with a different one
         * @param origID Identifier for the original target node
         * @param newID Identifier for the new target node
         */
        void replaceLeftArc(NodeID origID, NodeID newID) {
                if (nodeID > 0) {
                        if (dsNode->leftArc(origID) == NULL)
                                std::cout << "Paniek ! " << std::endl;
                } else
                        if (dsNode->rightArc(-origID) == NULL)
                                std::cout << "Paniek 2 ! " << std::endl;
                if (nodeID > 0) {
                        dsNode->leftArc(origID)->setNodeID(newID);
                } else {
                        dsNode->rightArc(-origID)->setNodeID(-newID);
                }
        }

        /**
         * Replace a right arc with a different one
         * @param origID Identifier for the original target node
         * @param newID Identifier for the new target node
         */
        void replaceRightArc(NodeID origID, NodeID newID) {
                if (nodeID > 0)
                        dsNode->rightArc(origID)->setNodeID(newID);
                else
                        dsNode->leftArc(-origID)->setNodeID(-newID);
        }

        /**
         * Get the sequence of this node
         * @return stl string containing the sequence
         */
        std::string getSequence() const {
                std::string seq = dsNode->getSequence();
                if (nodeID < 0)
                        Nucleotide::revCompl(seq);

                return seq;
        }

        /**
         * Get a subsequence of this node
         * @param offset Start offset
         * @param len Length of node
         * @return stl string containing the sequence
         */
        std::string substr(size_t offset, size_t len) const {
                if (nodeID > 0)
                        return dsNode->substr(offset, len);
                else
                        return Nucleotide::getRevCompl(dsNode->substr(length() - len - offset, len));
        }

        /**
         * Get a nucleotide at a specified position ('-' for out-of bounds)
         * @param pos Position in the sequence
         * @return Nucleotide at specified position
         */
        char getNucleotide(NodePosition pos) const {
                // check for out-of-bounds
                if (pos >= length())
                        return '-';
                if (nodeID < 0)
                        return Nucleotide::getComplement(dsNode->getNucleotide(length() - pos - 1));
                else
                        return dsNode->getNucleotide(pos);
        }

        /**
         * Set the sequence of this node
         * @param str String containing only 'A', 'C', 'G' and 'T'
         */
        void setSequence(const std::string& str) {
                if (nodeID > 0)
                        dsNode->setSequence(str);
                else
                        dsNode->setSequence(Nucleotide::getRevCompl(str));
        }

        /**
         * Get the left kmer of this node
         * @return The left kmer of this node
         */
        Kmer getLeftKmer() const {
                const std::string& seq = getSequence();
                return Kmer(seq);
        }

        /**
         * Get the right kmer of this node
         * @return The right kmer of this node
         */
        Kmer getRightKmer() const {
                const std::string& seq = getSequence();
                return Kmer(seq, seq.size() - Kmer::getK());
        }

        /**
         * Get the leftmost nucleotide of this node
         * @return The leftmost nucleotide
         */
        char peekNucleotideLeft() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideLeft();
                return Nucleotide::getComplement(dsNode->peekNucleotideRight());
        }

        /**
         * Get the rightmost nucleotide of this node
         * @return The rightmost nucleotide
         */
        char peekNucleotideRight() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideRight();
                return Nucleotide::getComplement(dsNode->peekNucleotideLeft());
        }

        /**
         * Get the nucleotide at position k - 1
         * @return The nucleotide at position k - 1
         */
        char peekNucleotideMarginalLeft() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideMarginalLeft();
                return Nucleotide::getComplement(dsNode->peekNucleotideMarginalRight());
        }

        /**
         * Get the nucleotide at position size - k
         * @return The nucleotide at position size - k
         */
        char peekNucleotideMarginalRight() const {
                if (nodeID > 0)
                        return dsNode->peekNucleotideMarginalRight();
                return Nucleotide::getComplement(dsNode->peekNucleotideMarginalLeft());
        }

        void inheritRightArcs(SSNode& target) {
                // make sure the node has currently no right arcs
                assert(numRightArcs() == 0);

                // update the arc information for the connected nodes
                for (ArcIt it = target.rightBegin(); it != target.rightEnd(); it++) {
                        SSNode rightNode(it);
                        rightNode.replaceLeftArc(target.getNodeID(), nodeID);
                }

                // copy the arcs
                copyRightArcs(target);
        }

        /**
         * Set the static node pointer
         * @param nodes Pointer to the double stranded nodes
         */
        static void setNodePointer(DSNode *nodes_) {
                nodes = nodes_;
        }
};

// ============================================================================
// NODE CHAIN CLASS
// ============================================================================

class NodeChain : public std::vector<NodeID>
{
private:
        size_t count;   // number of times the nodechain was observed

public:
        /**
         * Default constructor
         */
        NodeChain() : count(0) {}

        /**
         * Constructor from an input vector
         * @param input Input vector
         */
        NodeChain(const std::vector<NodeID>& input) : count(1) {
                for (const auto& it : input)
                        push_back(it);
        }

        /**
         * Constructor from an input vector range [first, last[
         * @param first Iterator the first element
         * @param last Iterator past the last element
         */
        NodeChain(std::vector<NodeID>::const_iterator first,
                  std::vector<NodeID>::const_iterator last) : count(1)
        {
                if (last <= first)
                        return;
                reserve(last - first);
                for (auto it = first; it < last; it++)
                        push_back(*it);
        }

        /**
         * Increment the count by one
         */
        void incrCount() {
                count++;
        }

        /**
         * Get the count
         * @return count
         */
        size_t getCount() const {
                return count;
        }

        /**
         * Set the count
         * @param target Target count
         */
        void setCount(size_t target) {
                count = target;
        }

        /**
         * Operator < overloading
         * @param rhs Right hand side
         */
      /*  bool operator<(const NodeChain& rhs) const {
                if (empty())

                for (size_t i = 0; i < size(); i++)
                        if ((*this)[i] != rhs[i])
                                return (*this)[i] < rhs[i];
                return false;
        }*/

        /**
         * Operator < overloading
         * @param rhs Right hand side
         */
        bool operator==(const NodeChain& rhs) const {
                if (size() != rhs.size())
                        return false;
                for (size_t i = 0; i < size(); i++)
                        if ((*this)[i] != rhs[i])
                                return false;
                return true;
        }

        /**
         * Get the reverse complement of the node chain
         * @return The reverse complementary chain
         */
        NodeChain getReverseComplement() const {
                NodeChain copy = *this;

                std::reverse(copy.begin(), copy.end());
                for (size_t i = 0; i < copy.size(); i++)
                        copy[i] = -copy[i];

                return copy;
        }

        /**
         * Get the representative node chain
         * @return The representative node chain
         */
        NodeChain getRepresentative() const {
                NodeChain RC = getReverseComplement();
                return (RC < *this) ? RC : *this;
        }

        void revCompl() {
                *this = getReverseComplement();
        }

        bool contains(const std::set<NodeRep>& noi) const {
                for (auto it : *this)
                        if (noi.find(NodeRep(it)) != noi.end())
                                return true;
                return false;
        }

        bool contains(NodeRep nr) const {
                for (auto it : *this)
                        if (NodeRep(it) == nr)
                                return true;
                return false;
        }

        bool contains(NodeID nodeID) const {
                for (auto it : *this)
                        if (it == nodeID)
                                return true;
                return false;
        }

        /**
         * Operator<< overloading
         * @param out Output file stream (input)
         * @param ncc NodeChainContainer to display
         * @return Output file stream
         */
        friend std::ostream &operator<<(std::ostream &out, const NodeChain& nc);

        static bool compCov(const NodeChain& l, const NodeChain& r) {
                return (l.getCount() < r.getCount());
        }

        /**
         * Sort by length (long to short), then lexigraphically
         */
        static bool compLenLex(const NodeChain& l, const NodeChain& r) {
                if (l.size() != r.size())
                        return (l.size() > r.size());
                for (size_t i = 0; i < l.size(); i++)
                        if (l[i] != r[i])
                                l[i] < r[i];
                return false;
        }

};

#endif
