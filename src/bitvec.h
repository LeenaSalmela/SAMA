/******************************************************************************
 *   Copyright (C) 2019 - 2020 Jan Fostier (jan.fostier@ugent.be)             *
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

#ifndef BITVEC_H
#define BITVEC_H

/**
 * The implementation implements the rank9 algorithm as described in
 * S. Vigna, "Broadword Implementation of Rank/Select Queries", WEA 2008
 * It relies on GCC's __builtin_popcountll, so please build this software
 * using the -mpopcnt flag to enable the SSE 4.2 POPCNT instruction.
 */

#include <vector>
#include <cstdint>

// ============================================================================
// BIT REFERENCE CLASS
// ============================================================================

class Bitref {

private:
        size_t& wordRef;        // reference to a word in the bitvector
        size_t bitmask;         // bitmask of the form (1 << bitIdx)

public:
        /**
         * Constructor
         * @param wordRef Reference to a word in the bitvector
         * @param bitmask Bitmask of the form (1 << bitIdx)
         */
        Bitref(size_t& wordRef, size_t bitmask) :
                wordRef(wordRef), bitmask(bitmask) {}

        /**
         * Set bitref to particular value
         * @param val Target value
         * @return Bitref reference after modification
         */
        const Bitref& operator=(bool val) {
                if (val)
                        wordRef |= bitmask;
                else
                        wordRef &= ~bitmask;
                return *this;
        }

        /**
         * Set bitref to another bitref
         * @param br Another bitref
         * @return Bitref reference
         */
        const Bitref& operator=(const Bitref& br) {
                return this->operator=(bool(br));
        }

        /**
         * Bool conversion operator
         */
        operator bool() const {
                return ((wordRef & bitmask) != 0) ? true : false;
        }
};

// ============================================================================
// BIT VECTOR CLASS
// ============================================================================

class Bitvec {

private:
        std::vector<size_t> data;       // actual bitvector
        std::vector<size_t> counts;     // interleaved 1st and 2nd level counts

public:
        /**
         * Get a bit at a certain position
         * @param p Position
         * @return true or false
         */
        bool operator[](size_t p) const {
                size_t w = p / 64;
                size_t b = p % 64;
                return (data[w] & (1ull << b)) != 0 ? true : false;
        }

        /**
         * Get a bit reference at a certain position
         * @param p Position
         * @return Bit reference object
         */
        Bitref operator[](size_t p) {
                size_t w = p / 64;
                size_t b = p % 64;
                return Bitref(data[w], 1ull << b);
        }

        /**
         * Create an index for the bitvector to support fast rank operations
         */
        void index() {
                counts = std::vector<size_t>((data.size() + 7) / 4, 0);

                size_t countL1 = 0, countL2 = 0;
                for (size_t w = 0, q = 0; w < data.size(); w++) {
                        if (w % 8 == 0) {       // store the L1 counts
                                countL1 += countL2;
                                counts[q] = countL1;
                                countL2 = __builtin_popcountll(data[w]);
                                q += 2;
                        } else {                // store the L2 counts
                                counts[q-1] |= (countL2 << (((w % 8)-1) * 9));
                                countL2 += __builtin_popcountll(data[w]);
                        }
                }
        }

        /**
         * Get the number of 1-bits within the range [0...p[ (preceding pos p)
         * @param p Position
         */
        size_t rank(size_t p) const {
                size_t w = p / 64;              // word index
                size_t q = (w / 8) * 2;         // counts index

                // add the first-level counts
                size_t rv = counts[q];

                // add the second-level counts
                int64_t t = (w % 8) - 1;
                rv += counts[q+1] >> (t + (t >> 60 & 8)) * 9 & 0x1FF;

                // add the popcount in the final word
                return rv + __builtin_popcountll( (data[w] << 1) << (63-(p % 64)) );
        }

        /**
         * Default Constructor
         */
        Bitvec() {}

        /**
         * Constructor
         * @param N Number of bits in the bitvector
         */
        Bitvec(size_t N) : data((N + 63) / 64, 0) {}
};

#endif
