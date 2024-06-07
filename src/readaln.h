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

#ifndef READALN_H
#define READALN_H

#include "settings.h"
#include "dbgraph.h"
#include "alignment.h"
#include "graphaln.h"
#include <mutex>
#include <sstream>

class GraphAlnWriter;

// ============================================================================
// GRAPH ALIGNMENT
// ============================================================================

class GraphAln : public std::vector<NodeID>
{
private:
        size_t firstNodeBegin;  // begin position in the first node
        size_t lastNodeEnd;     // end position of the last node
        size_t readBegin;       // begin position in the read
        size_t readEnd;         // end position of the aligned
        int score;              // alignment score

        /**
         * Debug routine: print information about a specific node
         * @param idx Index of the node
         * @param dBG Const-reference to de Bruijn graph
         */
        void printExt(size_t idx, const DBGraph& dBG) const;

        /**
         * Get the aligned node length for a node
         * @param idx Index of the node
         * @param dBG Const-reference to de Bruijn graph
         * @return The aligned node length
         */
        size_t getAlnLen(size_t idx, const DBGraph& dBG) const {
                assert(idx < size());
                if (size() == 1)        // special case: only one node
                        return lastNodeEnd - firstNodeBegin;
                if (idx == 0)           // special case: first node
                        return dBG.getSSNode(at(idx)).length() - firstNodeBegin;
                if (idx == size()-1)    // special case: last node
                        return lastNodeEnd;
                // general case
                return dBG.getSSNode(at(idx)).getMarginalLength();
        }

public:
        /**
         * Default constructor (empty alignment)
         */
        GraphAln() : firstNodeBegin(0), lastNodeEnd(0),
                readBegin(0), readEnd(0), score(0) {}

        /**
         * Construct a graph alignment on a single node
         * @param NodeID Node identifier
         * @param firstNodeBegin Begin position in the first node
         * @param lastNodeEnd End position of the last node
         * @param readBegin Begin position in the read
         * @param readEnd End position of the aligned
         * @param score Alignment score
         */
        GraphAln(NodeID nodeID, size_t firstNodeBegin, size_t lastNodeEnd,
                 size_t readBegin, size_t readEnd, int score) :
                 firstNodeBegin(firstNodeBegin), lastNodeEnd(lastNodeEnd),
                 readBegin(readBegin), readEnd(readEnd), score(score) {
                        push_back(nodeID);
        }

        /**
         * Extend the graph alignment by extending the present node
         * @param lastNodeEnd End position of the last node
         * @param readEnd End position of the aligned
         * @param score Alignment score
         */
        void extend(size_t lastNodeEnd, size_t readEnd, int score) {
                this->lastNodeEnd = lastNodeEnd;
                this->readEnd = readEnd;
                this->score = score;
        }

        /**
         * Extend the alignment by adding a new node
         * @param NodeID Node identifier
         * @param lastNodeEnd End position of the last node
         * @param readEnd End position of the aligned
         * @param score Alignment score
         */
        void extend(NodeID nodeID, size_t lastNodeEnd,
                    size_t readEnd, int score) {
                push_back(nodeID);
                extend(lastNodeEnd, readEnd, score);
        }

        /**
         * Extract the path sequence from the graph alignment
         * @param dBG Const-reference to de Bruijn graph
         * @return The path sequence from the graph alignment
         */
        std::string getPathStr(const DBGraph& dBG) {
                string result;

                for (size_t i = 0; i < size(); i++) {
                        SSNode n = dBG.getSSNode(at(i));

                        size_t b = (i == 0) ? firstNodeBegin : Kmer::getK() - 1;
                        size_t e = (i == size()-1) ? lastNodeEnd : n.length();

                        result.append(n.substr(b, e-b));
                }

                return result;
        }

        /**
         * Trim at most numChar characters from aln. Removes entire nodes.
         * Assumes that readEnd can be trimmed accordingly.
         * @param numChar Number of characters to remove from aln
         * @param dBG const-reference to de Bruijn graph
         * @param match Match score
         */
        void trimAln(size_t numChar, const DBGraph& dBG, int match) {
                size_t charTrimmed = 0;

                while (size() > 1) {    // leave at least one node
                        size_t nodeSize = lastNodeEnd + 1 - Kmer::getK();
                        if (charTrimmed + nodeSize > numChar)
                                return; // don't trim more than numChar

                        charTrimmed += nodeSize;
                        pop_back();     // remove final node
                        lastNodeEnd = dBG.getSSNode(back()).length();
                        readEnd -= nodeSize;
                        score -= nodeSize * match;
                }
        }

        /**
         * Reverse-complement the alignment
         * @param dBG Const-reference to de Bruijn graph
         * @param m Pattern length
         */
        void reverseComplement(const DBGraph& dBG, size_t m);

        /**
         * Print entire alignment
         * @param dBG Const-reference to de Bruijn graph
         */
        void print(const DBGraph& dBG) const {
                for (size_t i = 0; i < size(); i++)
                        printExt(i, dBG);
        }

        /**
         * Print only the last added node
         * @param dBG Const-reference to de Bruijn graph
         */
        void printLastNode(const DBGraph& dBG) const {
                if (empty())
                        return;
                printExt(size()-1, dBG);
        }

        /**
         * Return the begin position of the aligned read
         * @return The begin position of the aligned read
         */
        size_t getReadBegin() const {
                return readBegin;
        }

        /**
         * Get the end position of the read in the alignment
         * @return The end position of the read in the alignment
         */
        size_t getReadEnd() const {
                return readEnd;
        }

        /**
         * Get the length of the read in the alignment
         * @return The length of the read in the alignment
         */
        size_t getReadLen() const {
                return readEnd - readBegin;
        }

        /**
         * Get the begin position in the first node
         * @return The end begin position in the first node
         */
        size_t getFirstNodeBegin() const {
                return firstNodeBegin;
        }

        /**
         * Get the end position of the last node in the alignment
         * @return The end position of the last node in the alignment
         */
        size_t getLastNodeEnd() const {
                return lastNodeEnd;
        }

        /**
         * Get the alignment score
         * @return The alignment score
         */
        int getScore() const {
                return score;
        }

        /**
         * Set the alignment score
         * @param target The alignment score
         */
        void setScore(int target) {
                score = target;
        }

        /**
         * Get the alignment depth (= number of nodes in the alignment)
         * @return The alignment score
         */
        int getDepth() const {
                return (int)size();
        }

        /**
         * Write the alignment path to output file stream
         * @param ofs Opened output file stream
         */
        void write(std::ofstream& ofs) const {
                for (size_t i = 0; i < size(); i++) {
                        ofs << at(i);
                        if (i < size()-1)
                                ofs << "\t";
                }
                ofs << "\n";
        }

        /**
         * Write the alignment path to an open sequence file
         * @param file Opened output file stream
         */
        void write(SeqFile& file) const {
                std::stringstream ss;
                for (size_t i = 0; i < size(); i++) {
                        ss << at(i);
                        if (i < size()-1)
                                ss << "\t";
                }
                ss << "\n";
                file.writeLine(ss.str());
        }

        /**
         * Clear a graph alignment
         */
        void clear() {
                std::vector<NodeID>::clear();
                firstNodeBegin = lastNodeEnd = readBegin = readEnd = score = 0;
        }
};

// ============================================================================
// ALIGNMENT METRICS CLASS
// ============================================================================

class AlignmentMetrics
{
private:
        size_t numReads;                // number of reads handled
        size_t numCorrReads;            // number of reads corrected
        size_t numSubstitutions;        // number of substitutions
        size_t numIndels;               // number of indels
        size_t totReadLen;              // total read length
        size_t maxReadLen;              // maximum read length

        std::mutex metricMutex;         // mutex for merging metrics

public:
        /**
         * Default constructor
         */
        AlignmentMetrics() : numReads(0), numCorrReads(0), numSubstitutions(0),
                numIndels(0), totReadLen(0), maxReadLen(0) {}

        /**
         * Add the observation about a particular read
         * @param corrected True if (part of) the read was corrected
         * @param numSubstitutions Number of substitutions made to this read
         * @param numIndels Number of indels in this read
         * @param readLen Length of the read
         */
        void addObservation(bool corrected, size_t numSubstitutions,
                            size_t numIndels, size_t readLen) {

                if (corrected)
                        numCorrReads++;
                this->numSubstitutions += numSubstitutions;
                this->numIndels += numIndels;
                maxReadLen = max<size_t>(maxReadLen, readLen);
                totReadLen += readLen;
                numReads++;
        }

        /**
         * Merge metrics (thread-safe)
         * @param metrics Metrics to merge into this object
         */
        void mergeMetrics(const AlignmentMetrics& rhs);

        /**
         * Output statistics to the stdout
         */
        void printStatistics() const;
};

// ============================================================================
// READ CORRECTION CLASS
// ============================================================================

class ReadAligner
{
private:
        const DBGraph& dBG;             // const-reference to de Bruijn graph
        const Settings& settings;       // const-reference to settings object
        const KmerNPPTable& table;      // const-ref to <kmer, NodePos> table
        NWAligner alignment;            // Needleman-Wunsch aligner object
        NWAligner alignment2W;          // Needleman-Wunsch aligner object
        GraphAligner graphaln;          // Graph aligner object

        /**
         * Given a read, extract a set of seed (= partial alignments)
         * @param pattern Pattern to align
         * @param seeds Vector with consistent seeds (output)
         */
        void extractSeeds(const std::string& pattern,
                          std::vector<GraphAln>& seeds);

        /**
         * Find the node position pairs for a read
         * @param read Reference to the read
         * @param npp Vector of node position pairs
         */
        void findNPPFast(const std::string& read,
                         std::vector<NodePosPair>& npp);

        /**
         * Given a graph alignment, add all possible next nodes to the stack
         * @param P pattern to align
         * @param stack DFS stack <NodeID, depth> with next nodes to explore
         */
        void addAlnExt(const GraphAln& ga, const std::string& P,
                       std::vector<std::pair<NodeID, int> >& stack);

        /**
         * Extend a partial alignment using a DFS. The procedure extends the
         * alignment of P as far as adding a new nodes improves the alignment
         * score.
         * @param partAln Partial alignment of P, e.g. a seed alignment
         * @param P Pattern P to match
         * @return Search result (OPTIMAL, EXHAUSTED)
         */
        SearchRes extendAlnDFS(GraphAln& partAln, const std::string& P);

        /**
         * Extend and alignment both to the left and right
         * @param aln Partial alignment (seed) of pattern P
         * @param P Pattern P to align
         * @return Search result (OPTIMAL / EXHAUSTED)
         */
        SearchRes extendAln(GraphAln& aln, const std::string& P);

public:
        /**
         * Default constructor
         * @param dBG Reference to the De Bruijn graph
         * @param settings Reference to the settings class
         * @param table const-ref to <kmer, NodePos> table
         */
        ReadAligner(const DBGraph& dBG, const Settings& settings,
                    const KmerNPPTable& table) : dBG(dBG), settings(settings),
                    table(table), alignment(2, 1, -1, -3),
                    alignment2W(4, 1, -1, -3), graphaln(dBG) {}

        /**
         * Correct a specific read record
         * @param record Record to correct (input/output)
         * @param aln Resultin graph alignment (output)
         * @param metric Alignment metric to update (input/output)
         */
        void alignRead(FastQRecord& record, GraphAln& aln,
                       AlignmentMetrics& metric);
};

// ============================================================================
// READ ALIGNMENT HANDLER CLASS
// ============================================================================

class ReadAlnHandler
{
private:
        DBGraph& dBG;                   // const-ref to de Bruijn graph
        const Settings& settings;       // const-ref to settings object
        const KmerNPPTable& table;      // const-ref to <kmer, NodePos> table

        /**
         * Entry routine for worker thread
         * @param libaries Library container with libraries to be corrected
         * @param metrics Alignment metrics accross threads
         */
        void workerThread(FastQReader& myReader,
                          GraphAlnWriter& myWriter,
                          AlignmentMetrics& metrics);

public:
        /**
         * Default constructor
         * @param dBG const-reference to de Bruijn graph
         * @param settings const-reference to settings object
         * @param table const-ref to <kmer, NodePos> table
         */
        ReadAlnHandler(DBGraph& dBG, const Settings& settings,
                       const KmerNPPTable& table) : dBG(dBG),
                       settings(settings), table(table) {}

        /**
         * Align read libraries to de Bruijn graph
         * @param libraries Library container with reads
         */
        void align(LibraryContainer& libraries);
};

// ============================================================================
// GRAPH ALIGNMENT WRITER
// ============================================================================

class GraphAlnWriter
{
private:
        std::string filename1;          // name of the input file /1
        FileType fileType1;             // file type (FASTQ, FASTQ.GZ)

        std::string filename2;          // name of the input file /2
        FileType fileType2;             // file type (FASTQ, FASTQ.GZ)

        bool pairedEnd;                 // true for paired-end reads, false
                                        // for single-end reads (ignore file 2)

        std::thread oThread;            // output thread

        int maxNumChunks;

        // output thread variables (protect by outputMutex)
        std::mutex outputMutex;                 // output thread mutex
        std::condition_variable outputReady;    // output blocks ready condition
        std::condition_variable outputSpace;    // output blocks ready condition
        std::map<size_t, std::vector<GraphAln> > outputChunks;     // output blocks to be written
        size_t nextChunkID;                     // next chunk to be written

        /**
         * Entry routine for the output thread
         */
        void writerThread();

public:
        /**
         * Default constructor
         * @param filename1 name of the input file /1
         * @param filename2 name of the input file /2
         */
        GraphAlnWriter(const std::string& filename1, const std::string& filename2);

        /**
         * Start the writer thread
         * @param maxNumChunks Maximum number of chunks to buffer
         */
        void startWriterThread(size_t maxNumChunks);

        /**
         * Commit chunk of data to be written
         * @param chunk Buffer in which to store the records (output)
         * @param chunkID Unique chunk identifier (output)
         */
        void commitChunk(std::vector<GraphAln>& chunk, size_t chunkID);

        /**
         * Join the writer thread
         */
        void joinWriterThread();
};

#endif
