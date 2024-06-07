/******************************************************************************
 *   Copyright (C) 2014 - 2022 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of Detox    
 * 
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#ifndef FASTQ_H
#define FASTQ_H

#include "seqfile.h"

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <map>
#include <cassert>

#include <thread>
#include <mutex>
#include <condition_variable>

// ============================================================================
// FASTQ RECORD
// ============================================================================

class FastQRecord {

private:
        std::string seqID;      // sequence identifier
        std::string read;       // read itself
        std::string qual;       // quality string

public:
        /**
         * Clear the record
         */
        void clear() {
                seqID.clear();
                read.clear();
                qual.clear();
        }

        /**
         * Get a const-reference to the sequence identifier
         * @return a const-reference to the sequence identifier
         */
        const std::string& getSeqID() const {
                return seqID;
        }

        /**
         * Get a const-reference to the read
         * @return a const-reference to the read
         */
        const std::string& getRead() const {
                return read;
        }

        /**
         * Get a const-reference to the quality string
         * @return a const-reference to the quality string
         */
        const std::string& getQual() const {
                return qual;
        }

        /**
         * Read record from read file
         * @param Opened read file
         * @return True upon successful read, false otherwise
         */
        bool readFromFile(SeqFile& input);

        /**
         * Write record to file
         * @param Opened read file
         */
        void writeToFile(SeqFile& input) const;
};

// ============================================================================
// READ BLOCK CLASS
// ============================================================================

/**
 * A read block object contains a number of single-end or paired-end reads that
 * are read from file(s). In case of paired-end reads, the reads are stored
 * in an interleaved manner. In general, a readblock consists of multiple read
 * chunks.
 */

class ReadBlock : public std::vector<FastQRecord>
{
private:
        size_t nextChunkOffset; // offset of the next record to be read

public:
        /**
         * Construct read block from a vector of FastQ records (deep copy)
         * @param buffer A vector of FastQ records
         */
        ReadBlock(const std::vector<FastQRecord>& buffer) :
                std::vector<FastQRecord>(buffer), nextChunkOffset(0) {}

        /**
         * Default, deleted copy and default move constructor
         */
        ReadBlock() : nextChunkOffset(0) {}
        ReadBlock(const ReadBlock&) = delete;
        ReadBlock(ReadBlock&&) = default;

        /**
         * Deleted copy and default move assignement operator
         */
        ReadBlock& operator=(const ReadBlock&) = delete;
        ReadBlock& operator=(ReadBlock&& rhs) = default;


        /**
         * Return true of at least one chunk is available for reading
         * @return true or false
         */
        bool chunkAvailable() const {
                return nextChunkOffset < size();
        }

        /**
         * Get next available input record chunk
         * @param buffer Record buffer to write to (contents will be appended)
         * @param targetChunkSize Target chunk size (input)
         * @param pairedEnd If true, buffer.size() will always be even
         */
        void getNextChunk(std::vector<FastQRecord>& buffer,
                          size_t targetChunkSize, bool pairedEnd);

        /**
         * Read a block from file (paired-end reads)
         * @param file1 First fastq file
         * @param file2 First fastq file
         * @param targetBlockSize Desired number of nucleotides in this block
         */
        void readFromFile(SeqFile& file1, SeqFile& file2,
                          size_t targetBlockSize);

        /**
         * Read a block from file (paired-end reads)
         * @param file Fastq file
         * @param targetBlockSize Desired number of nucleotides in this block
         */
        void readFromFile(SeqFile& file1, size_t targetBlockSize);

        /**
         * Write a block to file (paired-end reads)
         * @param file1 First fastq file
         * @param file2 First fastq file
         */
        void writeToFile(SeqFile& file1, SeqFile& file2);

        /**
         * Write a block to file (paired-end reads)
         * @param file Fastq file
         */
        void writeToFile(SeqFile& file);
};

// ============================================================================
// FASTQ READER
// ============================================================================

class FastQReader
{
private:
        std::string filename1;          // name of the input file /1
        FileType fileType1;             // file type (FASTQ, FASTQ.GZ)
        std::string baseFilename1;      // file name w/o extension /1

        std::string filename2;          // name of the input file /2
        FileType fileType2;             // file type (FASTQ, FASTQ.GZ)
        std::string baseFilename2;      // file name w/o extension /1

        bool pairedEnd;                 // true for paired-end reads, false
                                        // for single-end reads (ignore file 2)

        const int numReadBlocks = 2;    // number of read blocks

        size_t targetChunkSize;         // target size for a single chunk
        size_t targetBlockSize;         // target size for a single block

        std::thread iThread;            // input thread

        // input thread variables (protect by inputMutex)
        std::mutex inputMutex;                  // input thread mutex
        std::condition_variable inputReady;     // input blocks ready condition
        std::vector<ReadBlock> inputBlocks;     // input blocks that are empty

        // worker thread variables (protect by workMutex)
        std::mutex workMutex;                   // worker thread mutex
        std::condition_variable workReady;      // work ready condition
        std::deque<ReadBlock> workBlocks;       // blocks being processed
        size_t chunkID;                         // unique ID per chunk

        /**
         * Entry routine for the input thread
         */
        void readerThread();

public:
        /**
         * Default constructor
         * @param filename1 name of the input file /1
         * @param filename2 name of the input file /2
         */
        FastQReader(const std::string& filename1, const std::string& filename2);

        /**
         * Return the filename(s) as a pair of strings
         * @return A pair of strings
         */
        std::pair<std::string, std::string> getFilename() const {
                return make_pair(filename1, filename2);
        }

        /**
         * Return the base filname as a pair of strings
         * @return A pair of strings
         */
        std::pair<std::string, std::string> getBaseFilename() const {
                return make_pair(baseFilename1, baseFilename2);
        }

        /**
         * Start the reader thread
         * @param targetChunkSize Target size for a single chunk
         * @param targetBlockSize Target size for a single block
         */
        void startReaderThread(size_t targetChunkSize, size_t targetBlockSize);

        /**
         * Get next read chunk from the input
         * @param chunk Buffer in which to store the records (output)
         * @param chunkID Unique chunk identifier (output)
         * @return False if no more reads are available, true otherwise
         */
        bool getNextChunk(std::vector<FastQRecord>& chunk, size_t& chunkID);

        /**
         * Join the reader thread
         */
        void joinReaderThread();
};

// ============================================================================
// FASTQ WRITER
// ============================================================================

class FastQWriter
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
        std::map<size_t, ReadBlock> outputChunks;     // output blocks to be written
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
        FastQWriter(const std::string& filename1, const std::string& filename2);

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
        void commitChunk(std::vector<FastQRecord>& chunk, size_t chunkID);

        /**
         * Join the writer thread
         */
        void joinWriterThread();
};

// ============================================================================
// LIBRARY CONTAINER CLASS
// ============================================================================

class LibraryContainer
{
private:
        typedef std::pair<std::string, std::string> pairOfStrings;

        // List of libraries stored as pairs-of-strings. In case of single-end
        // reads, the second filename is simply empty
        std::vector<pairOfStrings> filename;
        std::vector<pairOfStrings> baseFilename;

public:
        /**
         * Constructor
         * @param manifestFilename Read filename
         */
        LibraryContainer(const std::string& manifestFilename);

        /**
         * Get the number of libraries
         * @return The number of libraries
         */
        size_t size() const {
                return filename.size();
        }

        /**
         * Get the filename(s) for a specified library
         * @param index Input identifier
         * @return The filename(s) as a pair-of-strings
         */
        const pairOfStrings& getFilename(size_t index) const {
                return filename.at(index);
        }

        /**
         * Get the alignment filename(s) for a specified library
         * @param index Input identifier
         * @return The alignment filename(s) as a pair-of-strings
         */
        pairOfStrings getAlnFilename(size_t index) const {
                pairOfStrings base = baseFilename.at(index);
                return (base.second.empty()) ?
                        pairOfStrings(base.first + ".aln", std::string()) :
                        pairOfStrings(base.first + ".aln", base.second + ".aln");
        }

        /**
         * Get the base filename(s) for a specified library
         * @param index Input identifier
         * @return The base filename(s) as a pair-of-strings
         */
        const pairOfStrings& getBaseFilename(size_t index) const {
                return baseFilename.at(index);
        }
};

#endif
