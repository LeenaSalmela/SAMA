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

#include <functional>
#include "seqfile.h"
#include "fasta.h"

using namespace std;

// ============================================================================
// FASTA RECORD
// ============================================================================

bool FastARecord::readFromFile(SeqFile& inputFile)
{
        // read sequence identifier
        char c = inputFile.peekCharacter();

        // end of file might be reached
        if (!inputFile.good())
                return false;
        if (c != '>')
                throw ios::failure("File doesn't appear to be in FastA format");
        inputFile.getLine(seqID);
        if (!seqID.empty() && seqID.back() == '\n')
                seqID.pop_back();

        // read the actual read
        read.clear();
        while (inputFile.good() && inputFile.peekCharacter() != '>') {
                string line;
                inputFile.getLine(line);
                if (!line.empty() && line.back() == '\n')
                        line.pop_back();
                read.append(line);
        }

        return !read.empty();
}

void FastARecord::writeToFile(SeqFile& inputFile) const
{
        inputFile.writeLine(seqID);
        inputFile.writeChar('\n');
        inputFile.writeLine(read);
        inputFile.writeLine("\n");
}

// ============================================================================
// FASTA READER
// ============================================================================

FastAReader::FastAReader(const std::string& filename)
{
        // try to figure out the file format based on the ext
        string ext;
        FileType fileType;

        tie(fileType, std::ignore) = getFileType(filename);
        if ((fileType != FASTA) && (fileType != FASTA_GZ)) {
                string msg = "don't know how to open file: \"" + filename +
                             "\"\nExpected one of the following extensions: "
                             ".fasta or .fa (or .gz variants thereof)\n";
                throw runtime_error(msg);
        }

        if (!fileExists(filename))
                throw runtime_error("cannot open file " + filename);

        rf = move(SeqFile(fileType == FASTA_GZ));
        rf.open(filename);
}

bool FastAReader::getNextRecord(FastARecord& record) {
        return record.readFromFile(rf);
}
