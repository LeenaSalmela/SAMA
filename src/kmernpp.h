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

#ifndef KMERNPP_H
#define KMERNPP_H

#include "global.h"
#include "kmer/tkmer.h"
#include "util.h"

#include <utility>
#include <iostream>
#include <google/sparse_hash_map>

// ============================================================================
// PROTOTYPES
// ============================================================================

class DBGraph;

// ============================================================================
// NODE POSITION PAIR CLASS
// ============================================================================

class NodePosPair : public std::pair<NodeID, NodePosition>
{
public:
        /**
         * Default constructor
         */
        NodePosPair() : std::pair<NodeID, NodePosition>(0, 0) {}

        /**
         * Initialization constructor
         * @param id Node identifier
         * @param pos Position identifier
         */
        NodePosPair(NodeID id, NodePosition pos) :
                std::pair<NodeID, NodePosition>(id, pos) {}

        /**
         * Get the node identifier
         * @return Node identifier
         */
        NodeID getNodeID() const {
                return first;
        }

        /**
         * Reverse complement this node-position pair
         * @param nodeLength Length of the node ("marginal length" for k-mers,
         * "normal length" for actual positions)
         */
        NodePosPair getRevCompl(size_t nodeLength) const {
                return NodePosPair(-first, nodeLength - second - 1);
        }

        /**
         * Check whether the object points to valid position
         * @return True of false
         */
        bool isValid() const {
                return first != 0;
        }

        /**
         * Get the offset position
         * @return Offset position
         */
        NodePosition getPosition() const {
                return second;
        }

        /**
         * Set the offset position
         * @param target Target offset position
         */
        void setPosition(NodePosition target) {
                second = target;
        }

        /**
         * Operator << overloading
         * @param out Output stream to add NPP to
         * @param npp Right-hand side NPP
         * @return Output stream with the NPP added to it
         */
        friend std::ostream& operator<<(std::ostream &out, const NodePosPair& npp) {
                out << npp.first << " (" << npp.second << ")";
                return out;
        }
};

// ============================================================================
// <KMER - NODE-POSITION PAIR> TABLE
// ============================================================================

class KmerNPPTable {

private:
        const DBGraph& dBG;     // const-ref to de Bruijn graph
        bool doubleStranded;    // double stranded reads or not
        google::sparse_hash_map<Kmer, NodePosPair, KmerHash> table;

public:
        /**
         * Default constructor
         * @param dBG Const-reference to the Bruijn graph
         * @param doubleStranded Double stranded or not
         */
        KmerNPPTable(const DBGraph& dBG, bool doubleStranded) :
                dBG(dBG), doubleStranded(doubleStranded) {}

        /**
         * Resize hash table to desired size
         * @param size Desired size
         */
        void resize(size_t size) {
                table.resize(size);
        }

        /**
         * Clear the table
         */
        void clear() {
                table.clear();
        }

        /**
         * Insert a Kmer in the table
         * @param kmer Kmer to be inserted
         * @param npp Node-position pair
         * @return True if the kmer is inserted, false otherwise
         */
        bool insert(const Kmer &kmer, const NodePosPair& npp);

        /**
         * Write table to disk
         * @param filename Filename
         */
        void write(const std::string& filename) {
                FILE* myFile = fopen(filename.c_str(), "wb");
                table.write_metadata(myFile);
                table.write_nopointer_data(myFile);
                fclose(myFile);
        }

        /**
         * Load table from disk
         * @param filename Input filename
         * @return True if load was succesful, false otherwise
         */
        bool load(const std::string& filename) {
                if (!Util::fileExists(filename))
                        return false;

                FILE* myFile = fopen(filename.c_str(), "rb");
                table.read_metadata(myFile);
                table.read_nopointer_data(myFile);
                fclose(myFile);

                return true;
        }

        /**
         * Find a kmer in the table
         * @param kmer Kmer to look for
         * @return Node-position pair
         */
        NodePosPair find(const Kmer &kmer) const;
};

#endif
