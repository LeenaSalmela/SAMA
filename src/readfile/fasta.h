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

#ifndef FASTA_H
#define FASTA_H

#include "seqfile.h"
#include <string>

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class SeqFile;

// ============================================================================
// FASTA RECORD
// ============================================================================

class FastARecord {

private:
        std::string seqID;      // sequence identifier
        std::string read;       // read itself

public:
        /**
         * Clear the record
         */
        void clear() {
                seqID.clear();
                read.clear();
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
// FASTA READER
// ============================================================================

class FastAReader
{
private:
        std::string filename;           // name of the input file
        SeqFile rf;                     // read file

public:
        /**
         * Default constructor
         * @param filename name of the input file
         */
        FastAReader(const std::string& filename);

        /**
         * Return the filename
         * @return A string
         */
        std::string getFilename() const {
                return filename;
        }

        /**
         * Get next read chunk from the input
         * @param record FastA record (output)
         * @return False if no more reads are available, true otherwise
         */
        bool getNextRecord(FastARecord& record);
};

#endif
