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

#ifndef ARC_H
#define ARC_H

#include "global.h"

#include <atomic>
#include <fstream>
#include <cassert>
#include <vector>
#include <google/sparse_hash_map>

// ============================================================================
// EDGE
// ============================================================================

typedef std::pair<NodeID, NodeID> EdgeID;

// ============================================================================
// EDGE REPRESENTATIVE CLASS
// ============================================================================

/**
 * Represent an edge in a graph as a tuple (srcID, dstID). This class provides
 * convenient operator< and operator== overloads that take into account that
 * the edge (srcID, dstID) represents the same edge as (-dstID, -srcID)
 */
class EdgeRep {

private:
        NodeID srcID;           // source node identifier
        NodeID dstID;           // destination node identifier

public:
        /**
         * Constructor
         * @param srcID_ source node identifier
         * @param dstID_ destination node identifier
         */
        EdgeRep(NodeID srcID_, NodeID dstID_) {
                if (srcID_ < -dstID_) {
                        srcID = srcID_;
                        dstID = dstID_;
                } else {
                        srcID = -dstID_;
                        dstID = -srcID_;
                }
        }

        /**
         * Constructor
         * @param edge Edge specified as a tuple of nodes
         */
        EdgeRep(const std::pair<NodeID, NodeID>& edge) :
                EdgeRep(edge.first, edge.second) {}

        /**
         * Defaut constructor
         */
        EdgeRep() : srcID(0), dstID(0) {}

        /**
         * Get the source node identifier
         * @return the source node identifier
         */
        NodeID getSrcID() const {
                return srcID;
        }

        /**
         * Get the source node identifier
         * @return the source node identifier
         */
        NodeID getDstID() const {
                return dstID;
        }

        /**
         * Operator == overloading
         * @return true or false
         */
        bool operator==(const EdgeRep& rhs) const {
                return ((srcID == rhs.srcID) && (dstID == rhs.dstID));
        }

        /**
         * Operator < overloading
         * @return true or false
         */
        bool operator<(const EdgeRep& rhs) const {
                if (srcID != rhs.srcID)
                        return srcID < rhs.srcID;
                return dstID < rhs.dstID;
        }
        
        operator std::pair<NodeID, NodeID>() const {
                return std::pair<NodeID, NodeID>(srcID, dstID);
        }

        /**
         * Compute a hash function
         * @return Hash value
         */
        size_t getHash() const {

                size_t w = (size_t(srcID) << 32) + size_t(dstID);
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

        /**
         * Operator << overloading
         * @param out Output stream to add edge representative to
         * @param er Edge representative
         * @return Output stream with the edge representative added to it
         */
        friend std::ostream& operator<<(std::ostream &out, const EdgeRep& er) {
                out << "(" << er.srcID << " -> " << er.dstID << ")";
                return out;
        }
};

// ============================================================================
// HASH FUNCTION
// ============================================================================

struct EdgeHash {
        size_t operator()(const EdgeRep &er) const {
                return er.getHash();
        }
};

// ============================================================================
// EDGE MAP
// ============================================================================

template<class T>
using EdgeMap = google::sparse_hash_map<EdgeRep, T, EdgeHash>;

// ============================================================================
// ARC CLASS
// ============================================================================

class Arc {

private:
        NodeID nodeID;                  // ID of node to which the arc points
        std::atomic<Coverage> cov;      // arc coverage

public:
        /**
         * Default constructor
         */
        Arc() : nodeID(0), cov(0.0f) {}

        /**
         * Copy constructor
         * @param rhs Right hand side
         */
        Arc(const Arc& rhs) : nodeID(rhs.nodeID),
                cov(rhs.cov.load()) {}

        /**
         * Assignment operator
         * @param rhs Right hand side
         */
        Arc& operator=(const Arc& rhs) {
                nodeID = rhs.nodeID;
                cov = rhs.cov.load();
                return *this;
        }

        /**
         * Set the target nodeID the arc is pointing to
         * @param targetNodeID The ID of the node the arc is pointing to
         */
        void setNodeID(NodeID targetNodeID) {
                nodeID = targetNodeID;
        }

        /**
         * Get the target nodeID and the side of the target node
         * @return targetNodeID The ID of the node the arc is pointing to
         */
        NodeID getNodeID() const {
                return nodeID;
        }

        /**
         * Set the coverage of the arc
         * @param coverage Coverage of the arc
         */
        void setCov (Coverage coverage) {
                cov = coverage;
        }

        /**
         * Get the arc coverage
         * @return The arc coverage
         */
        Coverage getCov() const {
                return cov;
        }

        /**
         * Atomically increment the coverage
         * @param rhs Right hand side (default = 1.0f)
         */
        void incCov(Coverage rhs = 1.0f) {
                auto current = cov.load();
                while (!cov.compare_exchange_weak(current, current + rhs));
        }

        /**
         * Delete arc (mark as invalid)
         */
        void deleteArc() {
                nodeID = 0;
        }

        /**
         * Check if arc is valid (== not deleted)
         * @return True of false
         */
        bool isValid() const {
                return nodeID != 0;
        }

        /**
         * Write a node to file
         * @param ofs Open output file stream
         */
        void write(std::ofstream& ofs) const {
                ofs.write((char*)&nodeID, sizeof(nodeID));
                ofs.write((char*)&cov, sizeof(cov));
        }

        /**
         * Load a node from file
         * @param ifs Open input file stream
         */
        void read(std::ifstream& ifs) {
                ifs.read((char*)&nodeID, sizeof(nodeID));
                ifs.read((char*)&cov, sizeof(cov));
        }
};

// ============================================================================
// ARC ITERATOR CLASS
// ============================================================================

class ArcIt {

private:
        Arc* arcPtr;    // pointer to the positive arc
        bool reversed;  // true if this is an arc from a RC node
        Arc arc;        // arc copy, with nodeID sign swapped if necessary

        void setArc() {
                arc = *arcPtr;
                if (reversed)
                        arc.setNodeID(-arc.getNodeID());
        }

public:
        /**
         * Default constructor
         * @param id Identifier of the arc
         */
        ArcIt(Arc *arcPtr, bool reversed) : arcPtr(arcPtr), reversed(reversed) {
                setArc();
        }

        /**
         * Overloading of == operator
         * @return true of false
         */
        bool operator==(const ArcIt& it) {
                // make sure we're not comparing apples to pears
                assert(reversed == it.reversed);
                return arcPtr == it.arcPtr;
        }

        /**
         * Overloading the != operator
         * @return true of false
         */
        bool operator!=(const ArcIt &it) {
                // make sure we're not comparing apples to pears
                assert(reversed == it.reversed);
                return arcPtr != it.arcPtr;
        }

        /**
         * Overloading of postfix ++ operator
         * @return Copy of the iterator before the ++ operation
         */
        ArcIt operator++(int notused) {
                ArcIt copy = *this;
                arcPtr++;
                setArc();
                return copy;
        }

        /**
         * Overloading of prefix ++ operator
         * @return Reference to the iterator after ++ operator
         */
        ArcIt& operator++() {
                arcPtr++;
                setArc();
                return *this;
        }

        /**
         * Deference operator
         * @return a reference to the arc
         */
        const Arc& operator*() const {
                return arc;
        }

        /**
         * Deference operator
         * @return a reference to the arc
         */
        const Arc* operator->() const {
                return &arc;
        }
};

#endif
