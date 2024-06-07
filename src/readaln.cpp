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

#include <thread>
#include <string>
#include <iomanip>
#include <utility>
#include <algorithm>

#include "readaln.h"

using namespace std;

// ============================================================================
// GRAPH ALIGNMENT
// ============================================================================

void GraphAln::reverseComplement(const DBGraph& dBG, size_t m)
{
        GraphAln& aln = *this;          // short-hand notation

        if (empty())                    // if empty do nothing
                return;

        // reverse node chain
        reverse(begin(), end());        // reverse order
        for (auto& nodeID : aln)        // flip sign nodeID
                nodeID = -nodeID;

        // recompute firstNodeBegin and lastNodeEnd
        std::swap<size_t>(firstNodeBegin, lastNodeEnd);
        firstNodeBegin = dBG.getSSNode(front()).length() - firstNodeBegin;
        lastNodeEnd = dBG.getSSNode(back()).length() - lastNodeEnd;

        size_t readLen = readEnd - readBegin;
        readBegin = m - readEnd;
        readEnd = readBegin + readLen;
}

void GraphAln::printExt(size_t idx, const DBGraph& dBG) const
{
        assert(idx < size());

        const GraphAln& aln = *this;    // short-hand notation

        SSNode node = dBG.getSSNode(aln[idx]);

        size_t nodeBegin = (idx == 0) ? firstNodeBegin : Kmer::getK() - 1;
        size_t nodeEnd = (idx == size()-1) ? lastNodeEnd : node.length();

        string seq = node.getSequence().substr(nodeBegin, nodeEnd - nodeBegin);

        size_t whiteSp = getReadBegin();
        for (size_t i = 0; i < idx; i++)
                whiteSp += getAlnLen(i, dBG);

        for (size_t i = 0; i < whiteSp; i++)
                cout << " ";

        cout << seq << " (Node: " << aln[idx]
             << " [" << nodeBegin << ", " << nodeEnd << "[, score: "
             << aln.getScore() << ")\n";
}

// ============================================================================
// ALIGNMENT METRICS CLASS
// ============================================================================

void AlignmentMetrics::mergeMetrics(const AlignmentMetrics& rhs)
{
        lock_guard<mutex> lock(metricMutex);

        numCorrReads += rhs.numCorrReads;
        numSubstitutions += rhs.numSubstitutions;
        numIndels += rhs.numIndels;
        maxReadLen = max<size_t>(maxReadLen, rhs.maxReadLen);
        totReadLen += rhs.totReadLen;
        numReads += rhs.numReads;
}

void AlignmentMetrics::printStatistics() const
{
        size_t numUncorrected = numReads - numCorrReads;

        cout << "Read alignment report:\n";
        cout << "\tNumber of reads handled: " << numReads << endl;
        cout << "\tNumber of aligned reads: " << numCorrReads
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numCorrReads, numReads) << "%)" << endl;
        cout << "\tNumber of unaligned reads: " << numUncorrected
             << fixed << setprecision(2) << " ("
             << Util::toPercentage(numUncorrected, numReads) << "%)" << endl;
        cout << "\tNumber of substitutions in reads: " << numSubstitutions
             << fixed << setprecision(2) << " (avg of "
             << double(numSubstitutions)/double(numReads) << " per read)" << endl;
        cout << "\tNumber of indels in reads: " << numIndels
             << fixed << setprecision(2) << " (avg of "
             << double(numIndels)/double(numReads) << " per read)" << endl;
}

// ============================================================================
// READ CORRECTION CLASS
// ============================================================================

void ReadAligner::extractSeeds(const std::string& pattern,
                               vector<GraphAln>& seeds)
{
        // find the node position pairs using the kmer lookup table
        vector<NodePosPair> nppv(pattern.length() + 1 - Kmer::getK());
        findNPPFast(pattern, nppv);

        vector<bool> gappedSeed;
        size_t ip = nppv.size();        // previous index i
        for (size_t i = 0; i < nppv.size(); i++) {
                // consider two adjacent npps
                NodePosPair curr = nppv[i];
                NodePosPair prev = (ip != nppv.size()) ?
                        nppv[ip] : NodePosPair(0, 0);

                if (!curr.isValid())    // skip invalid npp
                        continue;

                size_t nb = nppv[i].getPosition();
                size_t ne = nb + Kmer::getK();
                size_t rb = i;
                size_t re = i + Kmer::getK();

                if (dBG.consecutiveNPP(prev, curr, i-ip)) {
                        if (curr.getPosition() == (prev.getPosition() + i - ip))
                                seeds.back().extend(ne, re, 0);
                        else
                                seeds.back().extend(curr.getNodeID(),
                                                    ne, re, 0);
                        if (i-ip > 1)   // flag gapped seed
                                gappedSeed.back() = true;
                } else {
                        seeds.push_back(GraphAln(curr.getNodeID(),
                                                 nb, ne, rb, re, 0));
                        gappedSeed.push_back(false);
                }

                ip = i;         // save the previous valid npp
        }

        // compute the alignment scores of the seeds
        const int M = alignment.getMatchScore();
        for (size_t i = 0; i < seeds.size(); i++) {
                GraphAln& s = seeds[i];

                // for ungapped seeds (= perfect matches), skip alignments
                if (!gappedSeed[i]) {
                        s.setScore(s.getReadLen() * M);
                        continue;
                }

                // perform an explicit alignment to compute the score
                string path = s.getPathStr(dBG);
                string Psub = pattern.substr(s.getReadBegin(), s.getReadLen());

                int score = alignment.alignBanded(path, Psub);
                s.setScore(score);
        }
}

void ReadAligner::findNPPFast(const string& read, vector<NodePosPair>& nppv)
{
        for (KmerIt it(read); it.isValid(); it++) {
                Kmer kmer = it.getKmer();
                NodePosPair npp = table.find(kmer);;
                nppv[it.getOffset()] = npp;

                if (!npp.isValid())
                        continue;

                // try to extend the exact match along the node
                NodeID nodeID = npp.getNodeID();
                const SSNode node = dBG.getSSNode(nodeID);

                size_t readPos = it.getOffset() + Kmer::getK();
                size_t nodePos = npp.getPosition() + Kmer::getK();

                while ((readPos < read.size()) && (nodePos < node.length())) {
                        if (read[readPos] != node.getNucleotide(nodePos))
                                break;

                        it++;
                        size_t nodeStart = nodePos + 1 - Kmer::getK();
                        nppv[it.getOffset()] = NodePosPair(nodeID, nodeStart);
                        nodePos++;
                        readPos++;
                }
        }
}

void ReadAligner::addAlnExt(const GraphAln& ga, const string& P,
                            vector<pair<NodeID, int> >& stack)
{
        // cannot extend empty alignment
        if (ga.empty())
                return;

        SSNode node = dBG.getSSNode(ga.back());

        char c = P[ga.getReadEnd()];
        NodeID addLast = 0;

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                NodeID nextID = it->getNodeID();
                SSNode next = dBG.getSSNode(it->getNodeID());
                if (next.peekNucleotideMarginalLeft() == c)
                        addLast = it->getNodeID();
                else
                        stack.push_back(make_pair(nextID, ga.getDepth() + 1));
        }

        if (addLast != 0)       // push delayed node = most promising node
                stack.push_back(make_pair(addLast, ga.getDepth() + 1));
}

SearchRes ReadAligner::extendAlnDFS(GraphAln& aln, const string& P)
{
        const int initDepth = aln.getDepth();
        const int initReadEnd = aln.getReadEnd();
        string Psub = P.substr(initReadEnd);

        alignment.reserveBanded(Psub.size(), aln.getScore());

        GraphAln& bestAln = aln;        // keep track of the best alignment
        GraphAln currAln = aln;         // keep track of the current alignment

        vector<int> offsetY;            // vertical offset in alignment matrix
        vector<pair<NodeID, int> > stack;       // todo stack <NodeID, depth>

        // extend the alignment of the last node before commencing a DFS
        SSNode node = dBG.getSSNode(currAln.back());

        string nodeStr = node.getSequence().substr(currAln.getLastNodeEnd());
        AlnRes2 res = alignment.alignBandedContd(Psub, nodeStr, 0);

        currAln.extend(currAln.getLastNodeEnd() + res.lenY,
                       initReadEnd + res.lenX, res.score);
        offsetY.push_back(nodeStr.size());

        // we've found a better alignment
        if (currAln.getScore() > bestAln.getScore())
                bestAln = currAln;

        // only continue further if you can possible improve the aln
        if (res.maxAtt > bestAln.getScore())
                addAlnExt(currAln, P, stack);

        size_t counter = 0;
        while (!stack.empty()) {
                if (counter++ > settings.getMaxDFSNodeVisited())
                        return EXHAUSTED;

                NodeID currID = stack.back().first;
                int depth = stack.back().second;
                stack.pop_back();

                // backtrack to the parent depth
                currAln.resize(depth - 1);
                offsetY.resize(depth - initDepth);

                // apply current node to the alignment
                SSNode currNode = dBG.getSSNode(currID);
                string nodeStr = currNode.getSequence().substr(Kmer::getK()-1);
                AlnRes2 res = alignment.alignBandedContd(Psub, nodeStr, offsetY.back());
                currAln.extend(currID, Kmer::getK() - 1 + res.lenY - offsetY.back(),
                               initReadEnd + res.lenX, res.score);

                offsetY.push_back(offsetY.back() + nodeStr.length());
                //currAln.printLastNode(dBG);

                // we've found a better alignment
                if (currAln.getScore() > bestAln.getScore())
                        bestAln = currAln;

                // only continue further if you can possible improve the aln
                if (res.maxAtt > bestAln.getScore())
                        addAlnExt(currAln, P, stack);
        }

        return OPTIMAL;
}

SearchRes ReadAligner::extendAln(GraphAln& aln, const string& P)
{
        // remove at most k-1 nucleotides from the seed
        aln.trimAln(Kmer::getK()-1, dBG, alignment.getMatchScore());
        SearchRes rightRes = extendAlnDFS(aln, P);

        string Prc = Nucleotide::getRevCompl(P);
        aln.reverseComplement(dBG, Prc.size());

        // remove at most k-1 nucleotides from the seed
        aln.trimAln(Kmer::getK()-1, dBG, alignment.getMatchScore());
        SearchRes leftRes = extendAlnDFS(aln, Prc);

        aln.reverseComplement(dBG, Prc.size());
        return ((rightRes == OPTIMAL) && (leftRes == OPTIMAL)) ?
                OPTIMAL : EXHAUSTED;
}

void ReadAligner::alignRead(FastQRecord& record, GraphAln& bestAln,
                            AlignmentMetrics& metrics)
{
        bestAln.clear();
        const string& read = record.getRead();

        // convert the <Node, Position> pairs into seeds (= graph alignments)
        vector<GraphAln> seeds;
        extractSeeds(read, seeds);

        // extend seeds by finding the best banded NW alignment on the graph
        for (GraphAln& aln : seeds) {
                //aln.print(dBG);

                SearchRes res = extendAln(aln, read);

                // save the best alignment
                if (aln.getScore() > bestAln.getScore())
                        bestAln = aln;
        }

        string path = bestAln.getPathStr(dBG);
        string Paln = read.substr(bestAln.getReadBegin(), bestAln.getReadLen());

        // use an aligner which supports twice as many indels (left + right)
        size_t nMatch, nSubst, nIndel;
        int score = alignment2W.alignBanded(path, Paln);
        alignment2W.getAlnStats(path, Paln, nMatch, nSubst, nIndel);

        bool corrected = !bestAln.empty();
        metrics.addObservation(corrected, nSubst, nIndel, read.length());
}

// ============================================================================
// READ ALIGNMENT HANDLER CLASS
// ============================================================================

void ReadAlnHandler::workerThread(FastQReader& myReader,
                                  GraphAlnWriter& myWriter,
                                  AlignmentMetrics& metrics)
{
        ReadAligner ra(dBG, settings, table);

        // local storage of reads
        vector<FastQRecord> readBuf;
        vector<GraphAln> writeBuf;

        // performance counters per thread
        AlignmentMetrics threadMetrics;

        size_t chunkID;
        while (myReader.getNextChunk(readBuf, chunkID)) {
                writeBuf.resize(readBuf.size());
                for (size_t i = 0; i < readBuf.size(); i++)
                        ra.alignRead(readBuf[i], writeBuf[i], threadMetrics);
                myWriter.commitChunk(writeBuf, chunkID);

        }

        // update the global metrics with the thread info (thread-safe)
        metrics.mergeMetrics(threadMetrics);
}

void ReadAlnHandler::align(LibraryContainer& libraries)
{
        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        const size_t ws = settings.getThreadIOWorkSize();
        for (size_t i = 0; i < libraries.size(); i++) {
                string fn1, fn2;
                tie(fn1, fn2) = libraries.getFilename(i);

                FastQReader myReader(fn1, fn2);
                if (fn2.empty()) {
                        cout << "Aligning single-end reads in file "
                             << fn1 << " to de Bruijn graph..." << endl;
                } else {
                        cout << "Aligning paired-end reads in files \"" << fn1
                             << "\" and \"" << fn2 << "\" to de Bruijn graph..."
                             << endl;
                }

                string outfn1, outfn2;
                tie(outfn1, outfn2) = libraries.getAlnFilename(i);
                GraphAlnWriter myWriter(outfn1, outfn2);
                myReader.startReaderThread(ws, ws * settings.getNumThreads());
                myWriter.startWriterThread(settings.getNumThreads());
                AlignmentMetrics metrics;

                // start worker threads
                vector<thread> worker(numThreads);
                for (size_t i = 0; i < worker.size(); i++)
                        worker[i] = thread(&ReadAlnHandler::workerThread,
                                   this, ref(myReader), ref(myWriter), ref(metrics));

                // wait for worker threads to finish
                for_each(worker.begin(), worker.end(), mem_fn(&thread::join));

                if (outfn2.empty()) {
                        cout << "Wrote file \"" << outfn1 << "\"" << endl;
                } else {
                        cout << "Wrote files \"" << outfn1 << "\" and \""
                             << outfn2 << "\"" << endl;
                }

                myReader.joinReaderThread();
                myWriter.joinWriterThread();
                metrics.printStatistics();
        }
}

// ============================================================================
// GRAPH ALIGNMENT WRITER
// ============================================================================

GraphAlnWriter::GraphAlnWriter(const string& filename1, const string& filename2) :
        filename1(filename1), fileType1(UNKNOWN_FT), filename2(filename2),
        fileType2(UNKNOWN_FT), pairedEnd(!filename2.empty())
{

}

void GraphAlnWriter::writerThread()
{
        // open the read file(s)
        SeqFile file1(fileType1 == FASTQ_GZ);
        file1.open(filename1, WRITE);

        SeqFile file2(fileType2 == FASTQ_GZ);
        if (pairedEnd) file2.open(filename2, WRITE);

        const size_t& termMsg = numeric_limits<size_t>::max();

        while (true) {
                // A) wait until the next chunk is available
                unique_lock<mutex> outputLock(outputMutex);
                outputReady.wait(outputLock,
                        [this, termMsg]{return  !outputChunks.empty() &&
                                               ((outputChunks.begin()->first == nextChunkID) ||
                                                (outputChunks.begin()->first == termMsg));});

                // check for termination message
                if (outputChunks.begin()->first == termMsg)
                        break;

                // get the output chunk
                vector<GraphAln> chunk = move(outputChunks.begin()->second);
                outputChunks.erase(outputChunks.begin());
                nextChunkID++;
                outputSpace.notify_all();
                outputLock.unlock();

                // B) write the chunk (only this thread has access)
                // no thread is held at this point
                if (pairedEnd) {
                        for (size_t i = 0; i < chunk.size(); ) {
                                chunk[i++].write(file1);
                                chunk[i++].write(file2);
                        }
                } else {
                        for (const auto& ga : chunk)
                                ga.write(file1);
                }
        }

        file1.close();
        if (pairedEnd) file2.close();
}

void GraphAlnWriter::commitChunk(vector<GraphAln>& chunk, size_t chunkID)
{
        // wait until space becomes available to store the read block
        // there is *always* space to store the next-to-be-written block
        unique_lock<mutex> outputLock(outputMutex);
        outputSpace.wait(outputLock,
                [this, chunkID]{return (chunkID == nextChunkID) ||
                                       (outputChunks.size() < maxNumChunks-1);});

        // deep copy
        outputChunks[chunkID] = chunk;

        // signal the output thread if necessary
        if (chunkID == nextChunkID)
                outputReady.notify_one();
}

void GraphAlnWriter::startWriterThread(size_t maxNumChunks)
{
        nextChunkID = 0;
        this->maxNumChunks = maxNumChunks;

        // start writer thread
        oThread = thread(&GraphAlnWriter::writerThread, this);
}

void GraphAlnWriter::joinWriterThread()
{
        // send a termination message (empty block)
        unique_lock<mutex> outputLock(outputMutex);
        outputChunks[numeric_limits<size_t>::max()] = vector<GraphAln>();
        outputReady.notify_one();
        outputLock.unlock();

        oThread.join();
}
