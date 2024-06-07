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

#include "alignment.h"

#include <stdlib.h>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <limits>

using namespace std;

// ============================================================================
// ALIGNMENT CLASS
// ============================================================================

void NWAligner::reserveBanded(size_t l1, size_t l2)
{
        // check the dimensions of s1 and s2
        int thisMaxDim = max(l1, l2);
        int thisMinDim = min(l1, l2);
        if ((thisMaxDim - thisMinDim) > maxIndel) {
                ostringstream oss;
                oss << "Cannot align sequences with length " << l1 << " and "
                    << l2 << " as the maximum indel is " << maxIndel << "\n";
                throw runtime_error(oss.str());
        }

        size_t reqSize = (2*maxIndel+1) * (thisMaxDim+1);

        // reallocate memory if necessary
        if (reqSize > currSize) {
                currSize = reqSize;
                delete [] matrix;
                matrix = new int[currSize];
        }
}

void NWAligner::reserveBanded(size_t m, int startScore)
{
        NWAligner& F = *this;   // shorthand notation

        size_t reqSize = (2*maxIndel+1) * (m+maxIndel+1);

        // reallocate memory if necessary
        if (reqSize > currSize) {
                currSize = reqSize;
                delete [] matrix;
                matrix = new int[currSize];
        }

        // initialize the borders of the matrix
        for (int i = 0; i <= maxIndel; i++) {
                F(i, 0) = startScore + i * G;
                F(0, i) = startScore + i * G;
        }
}

int NWAligner::S(char a, char b) const
{
        if (a == b)
                return M;
        if ((a == 'N') || (b == 'N'))
                return M;
        return I;
}

AlnRes NWAligner::alnGlobFreeEndGap(const string& s1, const string& s2)
{
        NWAligner& F = *this;   // shorthand notation

        // perform a normal, end-to-end alignment
        alignBanded(s1, s2);

        // remove trailing gaps
        int i = s1.size();
        int j = s2.size();

        // remove gaps at the end of s2, if any
        while ((i > 0) && (j < i + maxIndel) && (F(i, j) == F(i-1, j) + G))
                i--;
        if (i < s1.size())
                return AlnRes(i, j, F(i, j));

        // remove gaps at the end of s1, if any
        while ((j > 0) && (j > i - maxIndel) && (F(i, j) == F(i, j-1) + G))
                j--;
        if (j < s2.size())
                return AlnRes(i, j, F(i, j));

        // there were no trailing gaps, return the global aln score
        return AlnRes(i, j, F(i, j));
}

AlnRes NWAligner::alnGlobFreeEndGapS2(const string& s1, const string& s2)
{
        NWAligner& F = *this;   // shorthand notation

        // quick exit in case of too many gaps  FIXME
        if (abs((int)s1.length() - (int)s2.length()) > maxIndel)
                return AlnRes(0, 0, G*(s1.length() + s2.length()));

        // perform a normal, end-to-end alignment
        alignBanded(s1, s2);

        // remove trailing gaps
        int i = s1.size();
        int j = s2.size();

        // remove gaps at the end of s2, if any
        while ((i > 0) && (j < i + maxIndel) && (F(i, j) == F(i-1, j) + G))
                i--;
        if (i < s1.size())
                return AlnRes(i, j, F(i, j));

        // there were no trailing gaps, return the global aln score
        return AlnRes(i, j, F(i, j));
}

AlnRes NWAligner::align(const string& s1, const string& s2)
{
        // reserve memory (+1 for first/row column with gap penalties)
        D.resize(s1.length() + 1, s2.length() + 1);

        // initialize the borders of the matrix
        for (int i = 0; i <= s1.length(); i++)
                D(i, 0) = i * G;

        // initialize the borders of the matrix
        for (int j = 0; j <= s2.length(); j++)
                D(0, j) = j * G;

        // initialize the rest of the bulk of the matrix
        for (int i = 1; i <= (int)s1.length(); i++) {
                for (int j = 1; j <= (int)s2.length(); j++) {
                        int diag = D(i-1, j-1) + S(s1[i-1], s2[j-1]);
                        int up   = D(i-1, j) + G;
                        int left = D(i, j-1) + G;

                        D(i, j) = max(diag, max(up, left));
                }
        }

        return AlnRes(s1.length(), s2.length(), D(s1.length(), s2.length()) );
}

AlnRes NWAligner::trimTrailingGaps(const AlnRes& alnRes) const
{
        // remove trailing gaps
        int i = alnRes.s1len;
        int j = alnRes.s2len;

        // remove gaps at the end of s2, if any
        while ( (i > 0) && (D(i, j) == D(i-1, j) + G) )
                i--;
        if (i < alnRes.s1len)
                return AlnRes(i, j, D(i, j));

        // remove gaps at the end of s1, if any
        while ( (j > 0) && (D(i, j) == D(i, j-1) + G) )
                j--;
        if (j < alnRes.s2len)
                return AlnRes(i, j, D(i, j));

        // there were no trailing gaps, return the global aln score
        return alnRes;
}

int NWAligner::alignBanded(const string& s1, const string& s2)
{
        NWAligner& F = *this;   // shorthand notation

        // allocate memory if necessary
        reserveBanded(s1.length(), s2.length());

        // initialize the borders of the matrix
        for (int i = 0; i <= maxIndel; i++) {
                F(i, 0) = i * G;
                F(0, i) = i * G;
        }

        // initialize the rest of the bulk of the matrix
        for (int i = 1; i <= (int)s1.length(); i++)
        {
                int startj = max(1, i - maxIndel);
                int endj = min((int)s2.length(), i + maxIndel);

                for (int j = startj; j <= endj; j++)
                {
                        int diag = F(i-1, j-1) + S(s1[i-1], s2[j-1]);
                        int up   = (j < i + maxIndel) ? F(i-1, j) + G : diag-1;
                        int left = (j > i - maxIndel) ? F(i, j-1) + G : diag-1;

                        F(i, j) = max(diag, max(up, left));
                }
        }

        return F(s1.length(), s2.length());
}

AlnRes2 NWAligner::alignBandedContd(const string& X, const string& Y, int offsetY)
{
        NWAligner& F = *this;   // shorthand notation

        // last row of the matrix that will be computed
        int iLast = min(offsetY + Y.size(), X.size() + maxIndel);

        // initialize the rest of the bulk of the matrix
        for (int i = offsetY + 1; i <= iLast; i++) {
                int jBegin = max(1, i - maxIndel);
                int jLast = min<int>(X.size(), i + maxIndel);

                for (int j = jBegin; j <= jLast; j++) {
                        int diag = F(i-1, j-1) + S(Y[i-1-offsetY], X[j-1]);
                        int up   = (j < i + maxIndel) ? F(i-1, j) + G : diag-1;
                        int left = (j > i - maxIndel) ? F(i, j-1) + G : diag-1;

                        F(i, j) = max(diag, max(up, left));
                }
        }

        // get the aln score (X fully aligned) at position (sRow, X.size())
        int score = numeric_limits<int>::min(), sRow = 0;

        int iB = max<int>(X.size() - maxIndel, offsetY + 1);
        for (int i = iB; i <= iLast; i++) {
                if (score < F(i, X.size())) {
                        score = F(i, X.size());
                        sRow = i;
                }
        }

        // get the partial aln score at position (iLast, sCol) + max. attainable
        int partScore = numeric_limits<int>::min(), sCol = 0;
        int maxAtt = numeric_limits<int>::min();

        int jBegin = max(1, iLast - maxIndel);
        int jLast = min<int>(X.size(), iLast + maxIndel);
        for (int j = jBegin; j <= jLast; j++) {
                if (partScore < F(iLast, j)) {
                        partScore = F(iLast, j);
                        sCol = j;
                }
                maxAtt = max<int>(maxAtt, F(iLast, j) + M*(X.size()-j));
        }

        /*cout << "Full alignment score: " << score
             << ", at coordinates (" << sRow << ", " << X.size() << ")\n";
        cout << "Partial score: " << partScore
             << ", at coordinates (" << iLast << ", " << sCol << ")\n";
        cout << "Attainable score: " << maxAtt << endl;*/

        // print the alignment matrix
        /*cout << "\t";
        for (int i = 0; i < X.size(); i++)
                cout << "\t" << X[i];
        cout << "\n";

        for (int i = 0; i <= iLast; i++) {
                int jBegin = max(0, i - maxIndel);
                int jEnd = min<int>(X.size() + 1, i + maxIndel + 1);

                if (i >= offsetY + 1)
                        cout << Y[i-1-offsetY];
                cout << "\t";

                for (int i = 0; i < jBegin; i++)
                        cout << "\t";
                for (int j = jBegin; j < jEnd; j++)
                        cout << F(i, j) << "\t";
                cout << "\n";
        }*/

        return (score > numeric_limits<int>::min()) ?
                AlnRes2(X.size(), sRow, score, maxAtt) :
                AlnRes2(sCol, iLast, partScore, maxAtt);
}

NWAligner::NWAligner(int maxIndel, int M, int I, int G) :
        currSize(100), maxIndel(maxIndel), M(M), I(I), G(G)
{
        matrix = new int[currSize];
}

void NWAligner::printMatrix(const AlnRes& alnRes) const
{
        for (int i = 0; i <= alnRes.s1len; i++) {
                for (int j = 0; j <= alnRes.s2len; j++)
                        cout << D(i, j) << "\t";
                cout << "\n";
        }

        /*for (int l = 0; l < 2*maxIndel+1; l++) {
                for (int k = 0; k < maxDim+1; k++)
                        cout << matrix[k * (2 * maxIndel + 1) + l] << "\t";
                cout << "\n";
        }*/
}

void NWAligner::printAlignment(const AlnRes& alnRes,
                               const string& s1, const string& s2) const
{
        string al1, al2, mid;

        int i = alnRes.s1len;
        int j = alnRes.s2len;
        while (i > 0 && j > 0) {
                if ((i > 0) && (D(i, j) == D(i-1, j) + G)) {
                        al1.push_back(s1[i-1]);
                        al2.push_back('-');
                        mid.push_back(' ');
                        i--;
                } else if ((j > 0) && (D(i, j) == D(i, j-1) + G)) {
                        al1.push_back('-');
                        al2.push_back(s2[j-1]);
                        mid.push_back(' ');
                        j--;
                } else {
                        al1.push_back(s1[i-1]);
                        al2.push_back(s2[j-1]);
                        char c = (s1[i-1] == s2[j-1]) ? '|' : '*';
                        mid.push_back(c);
                        i--;
                        j--;
                }
        }

        for (int c = 0; c < i; c++) {
                al2.push_back(' ');
                mid.push_back(' ');
        }

        for (int c = 0; c < j; c++) {
                al1.push_back(' ');
                mid.push_back(' ');
        }

        reverse(al1.begin(), al1.end());
        reverse(al2.begin(), al2.end());
        reverse(mid.begin(), mid.end());

        al1 = s1.substr(0, i) + al1 + s1.substr(alnRes.s1len);
        al2 = s2.substr(0, j) + al2 + s2.substr(alnRes.s2len);

        cout << "Overlap alignment X[" << i << "-" << alnRes.s1len-1 << "], Y["
             << j << "-" << alnRes.s2len-1 << "]\n";
        for (size_t i = 0; i < s1.size(); i += 80) {
                cout << al1.substr(i, 80) << "\n"
                     << mid.substr(i, 80) << "\n"
                     << al2.substr(i, 80) << "\n\n";
        }
        cout << "Alignment score: " << alnRes.score << endl;
}

void NWAligner::printAlignmentBanded(const string& s1, const string& s2) const
{
        const NWAligner& F = *this;   // shorthand notation

        string al1, al2;

        int i = s1.size();
        int j = s2.size();
        while (i > 0 || j > 0) {
                if ((i > 0) && (j > 0) && (F(i, j) == F(i-1, j-1) + S(s1[i-1], s2[j-1]))) {
                        al1.push_back(s1[i-1]);
                        al2.push_back(s2[j-1]);
                        i--;
                        j--;
                } else if ((i > 0) && (j < i + maxIndel) && (F(i, j) == F(i-1, j) + G)) {
                        al1.push_back(s1[i-1]);
                        al2.push_back('-');
                        i--;
                } else {
                        al1.push_back('-');
                        al2.push_back(s2[j-1]);
                        j--;
                }
        }

        reverse(al1.begin(), al1.end());
        reverse(al2.begin(), al2.end());

        cout << al1 << "\n";
        for (int i = 0; i < max((int)al1.size(), (int)al2.size()); i++)
                if (al1[i] == al2[i] || al1[i] == 'N' || al2[i] == 'N')
                        cout << "|";
                else
                        cout << "*";
        cout << "\n" << al2 << "\n";
}

void NWAligner::getAlnStats(const string& s1, const string& s2, size_t& nMatch,
                            size_t& nSubst, size_t& nIndel) const
{
        const NWAligner& F = *this;   // shorthand notation
        nMatch = nSubst = nIndel = 0;

        int i = s1.size();
        int j = s2.size();
        while (i > 0 || j > 0) {
                if ((i > 0) && (j > 0) && (F(i, j) == F(i-1, j-1) + S(s1[i-1], s2[j-1]))) {
                        if (s1[i-1] == s2[j-1])
                                nMatch++;
                        else
                                nSubst++;
                        i--;
                        j--;
                } else if ((i > 0) && (j < i + maxIndel) && (F(i, j) == F(i-1, j) + G)) {
                        nIndel++;
                        i--;
                } else {
                        nIndel++;
                        j--;
                }
        }
}

AlnRes NWAligner::overlapAln(const string &s1, const string &s2)
{
        // reserve memory (+1 for first/row column with gap penalties)
        D.resize(s1.length() + 1, s2.length() + 1);

        // initialize the borders of the matrix
        for (int i = 0; i <= s1.length(); i++)
                D(i, 0) = 0;

        // initialize the borders of the matrix
        for (int j = 0; j <= s2.length(); j++)
                D(0, j) = j * G;

        // initialize the rest of the bulk of the matrix
        int bestScore = numeric_limits<int>::min(), best_i, best_j;
        for (int i = 1; i <= (int)s1.length(); i++) {
                for (int j = 1; j <= (int)s2.length(); j++) {
                        int diag = D(i-1, j-1) + S(s1[i-1], s2[j-1]);
                        int up   = D(i-1, j) + G;
                        int left = D(i, j-1) + G;

                        D(i, j) = max(diag, max(up, left));

                        // find the best score in the final row/column
                        if ( (i == s1.length() || j == s2.length())
                                && (D(i, j) > bestScore) )
                        {
                                bestScore = D(i, j);
                                best_i = i;
                                best_j = j;
                        }
                }
        }

        return AlnRes(best_i, best_j, bestScore);
}
