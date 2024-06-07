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

#include <iostream>
#include <algorithm>
#include <functional>
#include <climits>

#include "fastq.h"

using namespace std;

// ============================================================================
// FASTQ RECORD
// ============================================================================

bool FastQRecord::readFromFile(SeqFile& inputFile)
{
        // read sequence identifier
        char c = inputFile.peekCharacter();

        // end of file might be reached
        if (!inputFile.good())
                return false;
        if (c != '@')
                throw ios::failure("File doesn't appear to be in FastQ format");
        inputFile.getLine(seqID);
        if (!seqID.empty() && seqID.back() == '\n')
                seqID.pop_back();

        // read the actual read
        inputFile.getLine(read);
        if (!read.empty() && read.back() == '\n')
                read.pop_back();

        // read the + line
        string dummy;
        inputFile.getLine(dummy);
        // read the quality scores
        inputFile.getLine(qual);
        if (!qual.empty() && qual.back() == '\n')
                qual.pop_back();

        return !read.empty();
}

void FastQRecord::writeToFile(SeqFile& inputFile) const
{
        inputFile.writeLine(seqID);
        inputFile.writeChar('\n');
        inputFile.writeLine(read);
        inputFile.writeLine("\n+\n");
        inputFile.writeLine(qual);
        inputFile.writeChar('\n');
}

// ============================================================================
// READ BLOCK CLASS
// ============================================================================

void ReadBlock::getNextChunk(vector<FastQRecord>& buffer,
                             size_t targetChunkSize, bool pairedEnd)
{
        buffer.clear();

        // find out how many records to copy
        size_t thisChunkSize = 0;
        while (nextChunkOffset < size()) {
                buffer.push_back(move((*this)[nextChunkOffset++]));
                thisChunkSize += buffer.back().getRead().size();
                if (pairedEnd) {
                        buffer.push_back(move((*this)[nextChunkOffset++]));
                        thisChunkSize += buffer.back().getRead().size();
                }
                if (thisChunkSize >= targetChunkSize)
                        break;
        }
}

void ReadBlock::readFromFile(SeqFile& file1, SeqFile& file2,
                             size_t targetBlockSize)
{
        // clear the block
        clear();
        nextChunkOffset = 0;

        // read in new contents
        FastQRecord recordA, recordB;
        size_t thisBlockSize = 0;

        while (thisBlockSize < targetBlockSize) {
                bool flag1 = recordA.readFromFile(file1);
                bool flag2 = recordB.readFromFile(file2);
                if (!flag1 || !flag2)
                        break;

                push_back(recordA);
                push_back(recordB);
                thisBlockSize += recordA.getRead().size();
                thisBlockSize += recordB.getRead().size();
        }

        if (file1.good() != file2.good())
                cerr << "WARNING: paired-end FastQ files contain different "
                     << "number of reads" << endl;
}

void ReadBlock::readFromFile(SeqFile& file, size_t targetBlockSize)
{
        // clear the block
        clear();
        nextChunkOffset = 0;

        // read in new contents
        FastQRecord record;
        size_t thisBlockSize = 0;

        while (thisBlockSize < targetBlockSize) {
                if (!record.readFromFile(file))
                        break;
                push_back(record);
                thisBlockSize += record.getRead().size();
        }
}

void ReadBlock::writeToFile(SeqFile& file1, SeqFile& file2)
{
        // make sure we have an even number of read records
        assert(size() % 2 == 0);

        for (size_t i = 0; i < size(); i += 2) {
                (*this)[ i ].writeToFile(file1);
                (*this)[i+1].writeToFile(file2);
        }
}

void ReadBlock::writeToFile(SeqFile& file)
{
        for (const auto& r : *this)
                r.writeToFile(file);
}

// ============================================================================
// FASTQ READER
// ============================================================================

FastQReader::FastQReader(const string& filename1, const string& filename2) :
        filename1(filename1), fileType1(UNKNOWN_FT), filename2(filename2),
        fileType2(UNKNOWN_FT), pairedEnd(!filename2.empty())
{
        // try to figure out the file format based on the ext
        string ext;

        tie(fileType1, baseFilename1) = getFileType(filename1);
        if (fileType1 == UNKNOWN_FT) {
                string msg = "don't know how to open file: \"" + filename1 +
                             "\"\nExpected one of the following extensions: "
                             ".fastq or .fq (or .gz variants thereof)\n";
                throw runtime_error(msg);
        }

        if (!fileExists(filename1))
                throw runtime_error("cannot open file " + filename1);

        // no second file specified, assume single-ended reads
        if (!pairedEnd)
                return;

        tie(fileType2, baseFilename2) = getFileType(filename2);
        if (fileType2 == UNKNOWN_FT) {
                string msg = "don't know how to open file: \"" + filename2 +
                             "\"\nExpected one of the following extensions: "
                             ".fastq or .fq (or .gz variants thereof)\n";
                throw runtime_error(msg);
        }

        if (!fileExists(filename2))
                throw runtime_error("cannot open file " + filename2);
}

void FastQReader::readerThread()
{
        // open the read file(s)
        SeqFile file1(fileType1 == FASTQ_GZ);
        file1.open(filename1);

        SeqFile file2(fileType2 == FASTQ_GZ);
        if (pairedEnd) file2.open(filename2);

        size_t numReads = 0;
        while (file1.good()) {
                // A) wait until an input block is available
                unique_lock<mutex> inputLock(inputMutex);
                inputReady.wait(inputLock, [this]{return !inputBlocks.empty();});

                // get the input block
                ReadBlock block = move(inputBlocks.back());
                inputBlocks.pop_back();
                inputLock.unlock();

                // B) fill up the block (only this thread has access)
                // no mutexes are held by this thread at this point
                if (pairedEnd)
                        block.readFromFile(file1, file2, targetBlockSize);
                else
                        block.readFromFile(file1, targetBlockSize);

                numReads += block.size();
                cout << "Processed " << numReads << " reads...\r";
                cout.flush();

                // empty block: end-of-file reached
                if (block.empty())
                        break;

                // C) push the record block onto the worker queue
                unique_lock<mutex> workLock(workMutex);
                workBlocks.push_back(move(block));

                // notify workers that more work is available
                workReady.notify_all();
        }

        file1.close();
        if (pairedEnd) file2.close();

        // send a termination message to the workers
        unique_lock<mutex> workLock(workMutex);
        workBlocks.push_back(ReadBlock());      // empty block == termination
        workReady.notify_all();

        cout << "Processed " << numReads << " reads   " << endl;
}

bool FastQReader::getNextChunk(vector<FastQRecord>& chunk, size_t& chunkID)
{
        // wait until work becomes available
        unique_lock<mutex> workLock(workMutex);
        workReady.wait(workLock, [this]{return !workBlocks.empty();});

        // check if it is a termination message (empty block)
        if (workBlocks.front().empty())
                return false;

        // get a chunk of work (contents of chunk will be overwritten)
        workBlocks.front().getNextChunk(chunk, targetChunkSize, pairedEnd);
        chunkID = this->chunkID++;

        // last chunk: move block from work queue back to the input queue
        if (!workBlocks.front().chunkAvailable()) {
                ReadBlock tmp = move(workBlocks.front());
                workBlocks.pop_front();
                workLock.unlock();

                unique_lock<mutex> inputLock(inputMutex);
                inputBlocks.push_back(move(tmp));
                inputReady.notify_one();
        }

        assert(!chunk.empty());
        return true;
}

void FastQReader::startReaderThread(size_t targetChunkSize,
                                    size_t targetBlockSize)
{
        this->targetChunkSize = targetChunkSize;
        this->targetBlockSize = targetBlockSize;
        chunkID = 0;

        // initialize (empty) blocks on the input stack
        inputBlocks.resize(numReadBlocks);

        // start reader thread
        iThread = thread(&FastQReader::readerThread, this);
}

void FastQReader::joinReaderThread()
{
        iThread.join();

        inputBlocks.clear();
        workBlocks.clear();
}

// ============================================================================
// FASTQ WRITER
// ============================================================================

FastQWriter::FastQWriter(const string& filename1, const string& filename2) :
        filename1(filename1), fileType1(UNKNOWN_FT), filename2(filename2),
        fileType2(UNKNOWN_FT), pairedEnd(!filename2.empty())
{
        // try to figure out the file format based on the ext
        string ext;

        tie(fileType1, ignore) = getFileType(filename1);
        if (fileType1 == UNKNOWN_FT) {
                string msg = "don't know how to write file: \"" + filename1 +
                             "\"\nExpected one of the following extensions: "
                             ".fastq or .fq (or .gz variants thereof)\n";
                throw runtime_error(msg);
        }

        // no second file specified, assume single-ended reads
        if (!pairedEnd)
                return;

        tie(fileType2, ignore) = getFileType(filename2);
        if (fileType2 == UNKNOWN_FT) {
                string msg = "don't know how to write file: \"" + filename2 +
                             "\"\nExpected one of the following extensions: "
                             ".fastq or .fq (or .gz variants thereof)\n";
                throw runtime_error(msg);
        }
}

void FastQWriter::writerThread()
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
                        [this, termMsg]{return !outputChunks.empty() &&
                                              ((outputChunks.begin()->first == nextChunkID) ||
                                               (outputChunks.begin()->first == termMsg));});

                // check for termination message
                if (outputChunks.begin()->first == termMsg)
                        break;

                // get the output chunk
                ReadBlock chunk = move(outputChunks.begin()->second);
                outputChunks.erase(outputChunks.begin());
                nextChunkID++;
                outputSpace.notify_all();
                outputLock.unlock();

                // B) write the chunk (only this thread has access)
                // no thread is held at this point
                if (pairedEnd)
                        chunk.writeToFile(file1, file2);
                else
                        chunk.writeToFile(file1);
        }

        file1.close();
        if (pairedEnd) file2.close();
}

void FastQWriter::commitChunk(vector<FastQRecord>& chunk, size_t chunkID)
{
        // wait until space becomes available to store the read block
        // there is *always* space to store the next-to-be-written block
        unique_lock<mutex> outputLock(outputMutex);
        outputSpace.wait(outputLock,
                [this, chunkID]{return (chunkID == nextChunkID) ||
                                       (outputChunks.size() < maxNumChunks-1);});

        // deep copy
        outputChunks[chunkID] = ReadBlock(chunk);

        // signal the output thread if necessary
        if (chunkID == nextChunkID)
                outputReady.notify_one();
}

void FastQWriter::startWriterThread(size_t maxNumChunks)
{
        nextChunkID = 0;
        this->maxNumChunks = maxNumChunks;

        // start writer thread
        oThread = thread(&FastQWriter::writerThread, this);
}

void FastQWriter::joinWriterThread()
{
        // send a termination message (empty block)
        unique_lock<mutex> outputLock(outputMutex);
        outputChunks[numeric_limits<size_t>::max()] = ReadBlock();
        outputReady.notify_one();
        outputLock.unlock();

        oThread.join();
}

// ============================================================================
// LIBRARY CONTAINER CLASS
// ============================================================================

LibraryContainer::LibraryContainer(const std::string& mf)
{
        try {                           // option A: mf is a FastQ file
                FastQReader reader(mf, "");
                filename.push_back(make_pair(mf, ""));
                baseFilename.push_back(reader.getBaseFilename());
        } catch (exception& e) {        // option B: mf is a manifest file
                ifstream ifs(mf);
                if (!ifs)
                        throw runtime_error("cannot open file " + mf);

                string line;
                while (getline(ifs, line)) {
                        if (line.empty())
                                continue;
                        const string tokens = " \t,;";
                        size_t first = line.find_first_of(tokens);
                        size_t last  = line.find_last_of(tokens);

                        // no delimiter -> single-ended reads
                        if (last >= line.size()) {
                                FastQReader reader(line, "");
                                filename.push_back(reader.getFilename());
                                baseFilename.push_back(reader.getBaseFilename());
                        } else {        // paired-end reads
                                string file1 = line.substr(0, first);
                                string file2 = line.substr(last+1);
                                FastQReader reader(file1, file2);
                                filename.push_back(reader.getFilename());
                                baseFilename.push_back(reader.getBaseFilename());
                        }
                }
        }
}
