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

#ifndef SEQFILE_H
#define SEQFILE_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef HAVE_ZLIB
        #include "zlib.h"
#endif

// ============================================================================
// ENUM TYPES
// ============================================================================

typedef enum { READ, WRITE } ReadFileMode;
typedef enum { FASTQ, FASTQ_GZ, FASTA, FASTA_GZ, UNKNOWN_FT } FileType;

std::ostream &operator<<(std::ostream &out, const FileType &fileType);

// ============================================================================
// AUXILIARY ROUTINES
// ============================================================================

std::pair<FileType, std::string> getFileType(const std::string& fn);
bool fileExists(const std::string& fn);

// ============================================================================
// READFILE HANDLER
// ============================================================================

class ReadFileHandler {

protected:
        ReadFileMode mode;

public:
        /**
         * Virtual destructor for your convenience
         */
        virtual ~ReadFileHandler() {};

        /**
         * Open a file with a given filename
         * @param filename File to open
         * @param mode Read (default) or write
         */
        virtual void open(const std::string& filename,
                          ReadFileMode mode = READ) = 0;

        /**
         * Check if the input is still valid
         * @return True of false
         */
        virtual bool good() = 0;

        /**
         * Read a line from file
         * @return Pointer to an internal buffer containing the line
         */
        virtual void getLine(std::string& line) = 0;

        /**
         * Write a line to file
         * @param line String to write
         */
        virtual void writeLine(const std::string &line) = 0;

        /**
         * Write a character to file
         * @param c Character to write
         */
        virtual void writeChar(char c) = 0;

        /**
         * Peek at the next character in the stream
         * @return The character
         */
        virtual char peekCharacter() = 0;

        /**
         * Get one character from the stream
         * @return The character
         */
        virtual char getCharacter() = 0;

        /**
         * Close the file
         */
        virtual void close() = 0;

        /**
         * Reset the inputfile
         */
        virtual void reset() = 0;
};

// ============================================================================
// REGULAR READFILE HANDLER
// ============================================================================

class RegularReadFileHandler : public ReadFileHandler {

protected:
        FILE *fh;                               // file handler
        const int bufSize = 1024;               // buffer size
        char *buffer;                           // buffer

public:
        RegularReadFileHandler() {
                buffer = (char*)malloc(bufSize * sizeof(char));
        }

        /**
         * Virtual destructor for you convenience
         */
        ~RegularReadFileHandler() {
                free(buffer);
        }

        /**
         * Open a file with a given filename
         * @param filename File to open
         * @param mode Read (default) or write
         */
        void open(const std::string& filename, ReadFileMode mode = READ);

        /**
         * Check if the input is still valid
         * @return True of false
         */
        bool good() {
                return !feof(fh);
        }

        /**
         * Read a line from file
         */
        void getLine(std::string& line) {
                line.clear();

                while (fgets(buffer, bufSize, fh) != NULL) {
                        line.append(buffer);
                        if (line.back() == '\n')
                                break;
                }
        }

        /**
         * Write a line to file
         * @param line String to write
         */
        void writeLine(const std::string &line) {
                fputs(line.c_str(), fh);
        }

        /**
         * Write a character to file
         * @param c Character to write
         */
        void writeChar(char c) {
               fputc(c, fh);
        }

        /**
         * Peek at the next charachter in the stream
         * @return The character
         */
        char peekCharacter() {
                char c = fgetc(fh);
                ungetc(c, fh);
                return c;
        }

        /**
         * Get one character from the stream
         * @return The character
         */
        char getCharacter() {
                char c = fgetc(fh);
                return c;
        }

        /**
         * Close the file
         */
        void close();

        /**
         * Reset the inputfile
         */
        void reset();
};

// ============================================================================
// GZIPPED READFILE HANDLER
// ============================================================================

#ifdef HAVE_ZLIB

class GZipReadFileHandler : public ReadFileHandler {

protected:
        gzFile ifs;                     // gzipped input file stream
        const int bufSize = 1024;       // buffer size
        char *buffer;                   // buffer

public:
        /**
         * Default constructor
         */
        GZipReadFileHandler() {
                buffer = (char*)malloc(bufSize * sizeof(char));
        }

        /**
         * Destructor
         */
        ~GZipReadFileHandler() {
                free(buffer);
        }

        /**
         * Open a file with a given filename
         * @param filename File to open
         * @param mode Read (default) of write
         * @return True upon succes, false otherwise
         */
        void open(const std::string& filename, ReadFileMode mode);

        /**
         * Check if the input is still valid
         * @return True of false
         */
        bool good() {
                return !gzeof(ifs);
        }

        void getLine(std::string& line)  {
                line.clear();

                while (gzgets(ifs, buffer, bufSize) != NULL) {
                        line.append(buffer);
                        if (line.back() == '\n')
                                break;
                }
        }

        /**
         * Write a line to file
         * @param line String to write
         */
        void writeLine(const std::string &line) {
                gzwrite(ifs, line.c_str(), line.length());
        }

        /**
         * Write a character to file
         * @param c Character to write
         */
        void writeChar(char c) {
                gzputc(ifs, c);
        }

        /**
         * Peek at the next charachter in the stream
         * @return The character
         */
        char peekCharacter() {
                char c;
                c = gzgetc(ifs);
                gzungetc(c, ifs);
                return c;
        }

        /**
         * Get one character from the stream
         * @return The character
         */
        char getCharacter() {
                return gzgetc(ifs);
        }

        /**
         * Close the file
         */
        void close();

        /**
         * Reset the inputfile
         */
        void reset();
};

#endif

// ============================================================================
// SEQUENCE FILE
// ============================================================================

class SeqFile {

protected:
        ReadFileHandler *rfHandler;     // read file handler

public:
        /**
         * Default constructor
         */
        SeqFile() : rfHandler(NULL) {};

        /**
         * Default constructor
         * @param gzipped True if the file is gzipped
         */
        SeqFile(bool gzipped);

        /**
         * Destructor
         */
        virtual ~SeqFile() {
                delete rfHandler;
        }

        /**
         * Default, deleted copy and default move constructor
         */
        SeqFile(const SeqFile&) = delete;
        SeqFile(SeqFile&& rhs) {
                std::swap(rfHandler, rhs.rfHandler);
        }

        /**
         * Deleted copy and default move assignement operator
         */
        SeqFile& operator=(const SeqFile&) = delete;
        SeqFile& operator=(SeqFile&& rhs) {
                std::swap(rfHandler, rhs.rfHandler);
                return *this;
        };

        /**
         * Open a file
         * @param filename File to open
         * @param mode READ or WRITE
         */
        void open(const std::string& filename, ReadFileMode mode = READ) {
                rfHandler->open(filename, mode);
        }

        /**
         * Read a line from file
         * @line Line that was read (empty if failed)
         */
        void getLine(std::string& line) {
                rfHandler->getLine(line);
        }

        /**
         * Write a line to file
         * @param line String to write
         */
        void writeLine(const std::string &line) {
                rfHandler->writeLine(line);
        }

        /**
         * Write a character to file
         * @param c Character to write
         */
        void writeChar(char c) {
                rfHandler->writeChar(c);
        }

        /**
         * Peek at the next character in the stream
         * @return The character
         */
        char peekCharacter() {
                return rfHandler->peekCharacter();
        }

        /**
         * Get one character from the stream
         * @return The character
         */
        char getCharacter() {
                return rfHandler->getCharacter();
        }

        /**
         * Check if the input is still valid
         * @return True of false
         */
        bool good() {
                return rfHandler->good();
        }

        /**
         * Reset file to starting position
         */
        void reset() {
                rfHandler->reset();
        }

        /**
         * Close a read file
         */
        void close() {
                rfHandler->close();
        }
};

#endif
