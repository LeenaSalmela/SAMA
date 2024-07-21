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

#ifndef DSNODE_H
#define DSNODE_H

#include "arc.h"
#include "kmer/tstring.h"

#include <map>
#include <set>
#include <atomic>

// ============================================================================
// DOUBLE STRANDED NODE CLASS
// ============================================================================

class DSNode {

private:
        static Arc* arcs;

        typedef union {
                struct Packed {
                        uint8_t numLeft:3;              // [0...4]
                        uint8_t numRight:3;             // [0...4]
                        uint8_t invalid:1;
                        uint8_t flag:1;
                } p;
                uint8_t up;
        } Bitfield;

        TString sequence;       // DNA sequence
        ArcID leftID;           // ID of the first left arc or merged node
        ArcID rightID;          // ID of the first right arc or merged node
        Bitfield arcInfo;       // number of arcs at each node

        std::atomic<Coverage> cov;
        std::atomic<bool> myFlag1;
        std::atomic<bool> myFlag2;
#ifndef AVG_COV
        std::atomic<Count> *counts;
        std::atomic<Count> *arcCounts;
#endif
public:
        /**
         * Set the static arc pointer
         * @param arcPtr The static arc pointer
         */
        static void setArcsPointer(Arc *arcPtr) {
                arcs = arcPtr;
        }

        /**
         * Default constructor
         */
        DSNode() : leftID(0), rightID(0), cov(0), myFlag1(false), myFlag2(false)
#ifndef AVG_COV
		 , counts(NULL), arcCounts(NULL)
#endif
        {
                arcInfo.up = 0;
        }
  
        /**
         * Destructor
         */
        ~DSNode() {
	  if (counts != NULL) {
	    delete [] counts;
	  }
	  if (arcCounts != NULL) {
	    delete [] arcCounts;
	  }
	}

        DSNode& operator=(DSNode&& rhs) {

                if (this == &rhs)
                        return *this;

                // copy/move all members
                sequence = std::move(rhs.sequence);
                leftID = rhs.leftID;
                rightID = rhs.rightID;
                arcInfo = rhs.arcInfo;
                cov = rhs.cov.load();
                myFlag1 = rhs.myFlag1.load();
                myFlag2 = rhs.myFlag2.load();
#ifndef AVG_COV
		initializeCounts();
		for(int i = 0; i < getMarginalLength(); i++) {
		  counts[i] = rhs.counts[i].load();
		}
		for(int i = 0; i < getMarginalLength()-1; i++) {
		  arcCounts[i] = rhs.arcCounts[i].load();
		}
#endif		
                // load rhs in an empty state
                rhs.leftID = rhs.rightID = 0;
                rhs.arcInfo.up = 0;
                rhs.cov = 0;
                rhs.myFlag1 = false;
                rhs.myFlag2 = false;
#ifndef AVG_COV
		delete [] rhs.counts;
		rhs.counts = NULL;
		delete [] rhs.arcCounts;
		rhs.arcCounts = NULL;
#endif
                return *this;
        }

#ifndef AVG_COV
        std::atomic<Count> *getCounts() {
	  return counts;
	}

        void clearCounts() {
	  if (counts != NULL) {
	    for(int i = 0; i < getMarginalLength(); i++) {
	      counts[i] = 0;
	    }
	  }
	}

        std::atomic<Count> *getArcCounts() {
	  return arcCounts;
	}

        void clearArcCounts() {
	  if (arcCounts != NULL) {
	    for(int i = 0; i < getMarginalLength()-1; i++) {
	      arcCounts[i] = 0;
	    }
	  }
	}

        /**
         * Set the count
	 * @param pos The position
         * @param target The count
         */
        void setCount(int pos, Count target) {
	  assert(pos < getMarginalLength() && pos >= 0);
	  counts[pos] = target;
        }

        /**
         * Get the count
	 * "param pos The position
         * @return The count
         */
        Count getCount(int pos) const {
	  assert(pos < getMarginalLength() && pos >= 0);
	  return counts[pos];
        }

        /**
         * Atomically increment the count
	 * @param pos The position
         * @param rhs Right hand side (default = 1.0f)
         */
        void incCount(int pos, Count rhs = 1) {
	  assert(pos < getMarginalLength() && pos >= 0);
	  auto current = counts[pos].load();
	  while (!counts[pos].compare_exchange_weak(current, current + rhs));
        }

        /**
         * Set the arc count
	 * @param pos The position
         * @param target The arc count
         */
        void setArcCount(int pos, Count target) {
	  assert(pos < getMarginalLength()-1 && pos >= 0);
	  arcCounts[pos] = target;
        }

        /**
         * Get the arc count
	 * "param pos The position
         * @return The arc count
         */
        Count getArcCount(int pos) const {
	  assert(pos < getMarginalLength()-1 && pos >= 0);
	  return arcCounts[pos];
        }

        /**
         * Atomically increment the arc count
	 * @param pos The position
         * @param rhs Right hand side (default = 1.0f)
         */
        void incArcCount(int pos, Count rhs = 1) {
	  assert(pos < getMarginalLength()-1 && pos >= 0);
	  auto current = arcCounts[pos].load();
	  while (!arcCounts[pos].compare_exchange_weak(current, current + rhs));
        }

#endif
        /**
         * Set the coverage
         * @param target The coverage
         */
        void setCov(Coverage target) {
                cov = target;
        }

        /**
         * Get the coverage
         * @return The coverage
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
         * Set the flag
         * @param flag True of false
         */
        void setFlag1(bool flag) {
                myFlag1 = flag;
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag1() const {
                return myFlag1;
        }

        /**
         * Set the flag
         * @param flag True of false
         */
        void setFlag2(bool flag) {
                myFlag2 = flag;
        }

        /**
         * Get the flag
         * @return True of false
         */
        bool getFlag2() const {
                return myFlag2;
        }

        void swapRightArcsSign() {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        arcs[rightID + i].setNodeID(-arcs[rightID + i].getNodeID());
        }

        void swapLeftArcsSign() {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        arcs[leftID + i].setNodeID(-arcs[leftID + i].getNodeID());
        }

        /**
         * Get the identifier for the first left arc
         * @return The identifier for the first left arc
         */
        ArcID getFirstLeftArcID() const {
                return leftID;
        }

        /**
         * Get the identifier for the first right arc
         * @return The identifier for the first right arc
         */
        ArcID getFirstRightArcID() const {
                return rightID;
        }

        /**
         * Set the identifier for the first left arc
         * @param The identifier for the first left arc
         */
        void setFirstLeftArcID(ArcID target) {
                leftID = target;
        }

        /**
         * Get the identifier for the first right arc
         * @param The identifier for the first right arc
         */
        void setFirstRightArcID(ArcID target) {
                rightID = target;
        }

        /**
         * Invalidate a node (= mark as deleted)
         */
        void invalidate() {
                arcInfo.p.invalid = 1;
                sequence.clear();
#ifndef AVG_COV
		if (counts != NULL) {
		  delete [] counts;
		  counts = NULL;
		}
		if (arcCounts != NULL) {
		  delete [] arcCounts;
		  arcCounts = NULL;
		}
#endif
        }

        /**
         * Check if a node is valid
         * @return True or false
         */
        bool isValid() const {
                return (arcInfo.p.invalid == 0);
        }

        /**
         * Get the length of the node (# nucleotides in DNA string)
         * @return The length of the node
         */
        size_t getLength() const {
                return sequence.getLength();
        }

        /**
         * The marginal length == length - k + 1
         * @return The marginal length of the node
         */
        size_t getMarginalLength() const {
                return getLength() - Kmer::getK() + 1;
        }

        /**
         * Set the number of left arcs
         * @param numLeft The number of left arcs
         */
        void setNumLeftArcs(uint8_t numLeft) {
                arcInfo.p.numLeft = numLeft;
        }

        /**
         * Set the number of right arcs
         * @param numright The number of right arcs
         */
        void setNumRightArcs(uint8_t numRight) {
                arcInfo.p.numRight = numRight;
        }

        /**
         * Get the number of left arcs
         * @return The number of left arcs
         */
        uint8_t numLeftArcs() const {
                return arcInfo.p.numLeft;
        }

        /**
         * Get the number of right arcs
         * @return The number of right arcs
         */
        uint8_t numRightArcs() const {
                return arcInfo.p.numRight;
        }

        /**
         * Delete all left arcs
         */
        void deleteLeftArcs() {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        arcs[leftID + i].deleteArc();
                arcInfo.p.numLeft = 0;
        }

        /**
         * Delete all right arcs
         */
        void deleteRightArcs() {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        arcs[rightID + i].deleteArc();
                arcInfo.p.numRight = 0;
        }

        /**
         * Get a specific left arc
         * @param nodeID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* leftArc(NodeID nodeID) {
                for (int i = 0; i < arcInfo.p.numLeft; i++)
                        if (arcs[leftID + i].getNodeID() == nodeID)
                                return arcs + leftID + i;
                return NULL;
        }

        /**
         * Get a specific right arc
         * @param nodeID Identifier for the target node
         * @return Pointer to the specific arc, NULL if not found
         */
        Arc* rightArc(NodeID nodeID) {
                for (int i = 0; i < arcInfo.p.numRight; i++)
                        if (arcs[rightID + i].getNodeID() == nodeID)
                                return arcs + rightID + i;
                return NULL;
        }

        /**
         * Delete a specific left arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteLeftArc(NodeID targetID);

        /**
         * Delete a specific right arc
         * @param targetID Identifier for the target node
         * @return True of the arc was deleted, false if the arc was not found
         */
        bool deleteRightArc(NodeID targetID);

        /**
         * Get an iterator pointing to the first left arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the first left arc
         */
        ArcIt leftBegin(bool reversed = false) const {
                return ArcIt(arcs + leftID, reversed);
        }

        /**
         * Get an iterator pointing past the last left arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the last left arc
         */
        ArcIt leftEnd(bool reversed = false) const {
                return ArcIt(arcs + leftID + arcInfo.p.numLeft, reversed);
        }

        /**
         * Get an iterator pointing to the first right arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the first left arc
         */
        ArcIt rightBegin(bool reversed = false) const {
                return ArcIt(arcs + rightID, reversed);
        }

        /**
         * Get an iterator pointing past the last right arc
         * @param reversed True if you want arcs from the negative node
         * @return An iterator pointing to the last right arc
         */
        ArcIt rightEnd(bool reversed = false) const {
                return ArcIt(arcs + rightID + arcInfo.p.numRight, reversed);
        }

        /**
         * Set the sequence of this node
         * @param str String containing only 'A', 'C', 'G' and 'T'
         */
        void setSequence(const std::string& str) {
                sequence.setSequence(str);
        }


#ifndef AVG_COV
        void initializeCounts() {
	  if (counts != NULL) {
	    delete [] counts;
	  }
	  counts = new std::atomic<Count>[getMarginalLength()];
	  for(int i = 0; i < getMarginalLength(); i++) {
	    counts[i] = 0;
	  }
	  
	  if (arcCounts != NULL) {
	    delete [] arcCounts;
	  }
	  arcCounts = new std::atomic<Count>[getMarginalLength()-1];
	  for(int i = 0; i < getMarginalLength()-1; i++) {
	    arcCounts[i] = 0;
	  }
	}
#endif

  
        /**
         * Get the sequence of this node
         * @return The sequence of this node
         */
        std::string getSequence() const {
                return sequence.getSequence();
        }

        /**
         * Get a subsequence of this node
         * @param offset Start offset
         * @param len Length of node
         * @return stl string containing the sequence
         */
        std::string substr(size_t offset, size_t len) const {
                return sequence.substr(offset, len);
        }

        /**
         * Get a nucleotide at a specified position ('-' for out-of bounds)
         * @param pos Position in the sequence
         * @return Nucleotide at specified position
         */
        char getNucleotide(NodePosition pos) const {
                // check for out-of-bounds
                if (pos >= getLength())
                        return '-';
                return sequence[pos];
        }

        /**
         * Get the tight sequence of this node
         * @return The tight sequence
         */
        const TString& getTSequence() const {
                return sequence;
        }

        /**
         * Get the leftmost nucleotide of this node
         * @return The leftmost nucleotide
         */
        char peekNucleotideLeft() const {
                return sequence.peekNucleotideLeft();
        }

        /**
         * Get the rightmost nucleotide of this node
         * @return The rightmost nucleotide
         */
        char peekNucleotideRight() const {
                return sequence.peekNucleotideRight();
        }

        /**
         * Get the nucleotide at position k - 1
         * @return The nucleotide at position k - 1
         */
        char peekNucleotideMarginalLeft() const {
                return sequence.peekNucleotideMarginalLeft();
        }

        /**
         * Get the nucleotide at position size - k
         * @return The nucleotide at position size - k
         */
        char peekNucleotideMarginalRight() const {
                return sequence.peekNucleotideMarginalRight();
        }

        /**
         * Get the leftmost kmer of this node
         * @return The leftmost kmer
         */
        Kmer getLeftKmer() const {
                return Kmer(sequence, 0);
        }

        /**
         * Get the rightmost kmer of this node
         * @return The rightmost kmer
         */
        Kmer getRightKmer() const {
                return Kmer(sequence, sequence.getLength() - Kmer::getK());
        }

        /**
         * Write a node to file
         * @param ofs Open output file stream
         */
        void write(std::ofstream& ofs) const {
                ofs.write((char*)&cov,sizeof(cov));
                ofs.write((char*)&leftID, sizeof(leftID));
                ofs.write((char*)&rightID, sizeof(rightID));
                ofs.write((char*)&arcInfo, sizeof(arcInfo));
                sequence.write(ofs);
#ifndef AVG_COV
		for(int i = 0; i < getMarginalLength(); i++) {
		  ofs.write((char*)&counts[i],sizeof(counts[i]));
		}
		for(int i = 0; i < getMarginalLength()-1; i++) {
		  ofs.write((char*)&arcCounts[i],sizeof(arcCounts[i]));
		}
#endif
        }

        /**
         * Load a node from file
         * @param ifs Open input file stream
         */
        void read(std::ifstream& ifs) {
                ifs.read((char*)&cov, sizeof(cov));
                ifs.read((char*)&leftID, sizeof(leftID));
                ifs.read((char*)&rightID, sizeof(rightID));
                ifs.read((char*)&arcInfo, sizeof(arcInfo));
                sequence.read(ifs);
#ifndef AVG_COV
		if (counts != NULL)
		  delete [] counts;
		counts = new std::atomic<Count>[getMarginalLength()];
		//std::cout << "ML: " << getMarginalLength() << std::endl;
		for(int i = 0; i < getMarginalLength(); i++) {
		  ifs.read((char*)&counts[i], sizeof(counts[i]));
		  //std::cout << "N: " <<  counts[i] << std::endl;
		}
		if (arcCounts != NULL)
		  delete [] arcCounts;
		arcCounts = new std::atomic<Count>[getMarginalLength()-1];
		for(int i = 0; i < getMarginalLength()-1; i++) {
		  ifs.read((char*)&arcCounts[i], sizeof(arcCounts[i]));
		  //std::cout << "A: " <<  counts[i] << std::endl;
		}
#endif
        }
};

#endif
