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

#include "dbgraph.h"
#include <iomanip>

using namespace std;

void DBGraph::sanityCheck()
{
        size_t countValidNodes = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode n = getSSNode(id);
                if (n.isValid())
                        countValidNodes++;
        }

        if (countValidNodes != numValidNodes)
                cerr << "Counted " << countValidNodes << " valid nodes, but "
                     << "internal variable has value: " << numValidNodes << endl;

        size_t countValidArcs = 0;
        for (ArcID id = 1; id <= numArcs; id++)
                if (arcs[id].isValid())
                        countValidArcs++;

        if (countValidArcs != numValidArcs)
                cerr << "Counted " << countValidArcs << " valid arcs, but "
                     << "internal variable has value: " << numValidArcs << endl;


        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)     // skip node 0, doesn't exist
                        continue;

                // shortcuts
                SSNode n = getSSNode(id);

                if (!n.isValid()) {      // invalid nodes shoud not have arcs
                        if ((n.numLeftArcs() > 0) || (n.numRightArcs() > 0))
                                cerr << "\t\tNode " << id
                                     << " is invalid but has arcs." << endl;
                        continue;
                }

                // check the continuity of connected nodes
                string sequence = n.getSequence();

                Kmer nLeft(sequence);
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        SSNode l = getSSNode(it->getNodeID());
                        if (!l.isValid())
                                cerr << "\t\tNode " << id << " has a left arc "
                                        "to an invalid node\n";

                        //assert(l.rightArc(id)->getCov() == it->getCov());
                        //assert(fabs(it->getCov()-l.rightArc(id)->getCov()) <= 1e-5);
                        if (l.rightArc(id)->getCov() != it->getCov())
                            cout << std::setprecision(8) << l.rightArc(id)->getCov() << " != " << it->getCov() << endl;

                        Kmer lRight = l.getRightKmer();
                        lRight.pushNucleotideRight(nLeft.peekNucleotideRight());
                        if (nLeft != lRight)
                                cerr << "\tError in contiguity of nodes\n";
                        assert(nLeft == lRight);
                        countValidArcs++;
                }

                Kmer nRight(sequence, sequence.size()-Kmer::getK());
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        SSNode r = getSSNode(it->getNodeID());
                        if (!r.isValid())
                                cerr << "\t\tNode " << id << " has a right arc "
                                        "to an invalid node\n";

                        //assert(r.leftArc(id)->getCov() == it->getCov());
                        //assert(fabs(it->getCov()-r.leftArc(id)->getCov()) <= 1e-5);
                        if (r.leftArc(id)->getCov() != it->getCov())
                            cout << std::setprecision(8) << r.leftArc(id)->getCov() << " != " << it->getCov() << endl;

                        Kmer rLeft = r.getLeftKmer();
                        rLeft.pushNucleotideLeft(nRight.peekNucleotideLeft());
                        if (nRight != rLeft)
                                cerr << "\tError in contiguity of nodes\n";
                        assert (nRight == rLeft);
                        countValidArcs++;
                }
        }
}
