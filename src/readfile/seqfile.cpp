/******************************************************************************
 *   Copyright (C) 2014 - 2022 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of Detox    
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

#include <cstring>
#include <cassert>
#include <algorithm>
#include "seqfile.h"

using namespace std;

// ============================================================================
// FILETYPE ENUM
// ============================================================================

std::ostream &operator<<(std::ostream &out, const FileType &fileType) {
        switch (fileType) {
                case FASTQ :
                        out << "fastq";
                        break;
                case FASTQ_GZ :
                        out << "fastq.gz";
                        break;
                case FASTA:
                        out << "fasta";
                        break;
                case FASTA_GZ:
                        out << "fasta.gz";
                        break;
                case UNKNOWN_FT:
                        out << "unknown";
                        break;
        }

        return out;
}

// ============================================================================
// AUXILIARY ROUTINES
// ============================================================================

pair<FileType, string> getFileType(const string& fn)
{
        string ext;

        if (fn.length() >= 3)
                ext = fn.substr(fn.length() - 3);
        transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
        if (ext == ".FQ")
                return make_pair(FASTQ, fn.substr(0, fn.length() - 3));
        else if (ext == ".FA")
                return make_pair(FASTA, fn.substr(0, fn.length() - 3));

        if (fn.length() >= 6)
                ext = fn.substr(fn.length() - 6);
        transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
        if (ext == ".FASTQ")
                return make_pair(FASTQ, fn.substr(0, fn.length() - 6));
        else if (ext == ".FQ.GZ")
                return make_pair(FASTQ_GZ, fn.substr(0, fn.length() - 6));
        else if (ext == ".FASTA")
                return make_pair(FASTA, fn.substr(0, fn.length() - 6));
        else if (ext == ".FA.GZ")
                return make_pair(FASTA_GZ, fn.substr(0, fn.length() - 6));

        if (fn.length() >= 9)
                ext = fn.substr(fn.length() - 9);
        transform(ext.begin(), ext.end(), ext.begin(), ::toupper);
        if (ext == ".FASTQ.GZ")
                return make_pair(FASTQ_GZ, fn.substr(0, fn.length() - 9));
        if (ext == ".FASTA.GZ")
                return make_pair(FASTA_GZ, fn.substr(0, fn.length() - 9));

        return make_pair(UNKNOWN_FT, string());
}

bool fileExists(const std::string& fn) {
        ifstream ifs(fn);
        return ifs.good();
}

// ============================================================================
// REGULAR READFILE HANLDER
// ============================================================================

void RegularReadFileHandler::open(const string& filename, ReadFileMode mode)
{
        this->mode = mode;
        if (mode == READ)
                fh = fopen(filename.c_str(), "r");
        else
                fh = fopen(filename.c_str(), "w");
}

void RegularReadFileHandler::reset()
{
        if (mode == READ)
                fseek(fh, 0, SEEK_SET);
}

void RegularReadFileHandler::close()
{
        fclose(fh);
}

// ============================================================================
// GZIPPED READFILE HANDLER
// ============================================================================

#ifdef HAVE_ZLIB

void GZipReadFileHandler::open(const std::string& filename, ReadFileMode mode)
{
        this->mode = mode;
        const char* m = (mode == READ) ? "r" : "w";
        ifs = gzopen(filename.c_str(), m);
}

void GZipReadFileHandler::reset() {
        gzrewind(ifs);
}

void GZipReadFileHandler::close()
{
        gzclose(ifs);
        ifs = Z_NULL;
}

#endif

// ============================================================================
// SEQUENCE FILE
// ============================================================================

SeqFile::SeqFile(bool gzipped) : rfHandler(NULL)
{
#ifndef HAVE_ZLIB
        if (gzipped)
                throw ios_base::failure("Software was compiled without "
                                        "zlib support");
        rfHandler = new RegularReadFileHandler();
#else
        if (gzipped)
                rfHandler = new GZipReadFileHandler();
        else
                rfHandler = new RegularReadFileHandler();
#endif
}
