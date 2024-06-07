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

#include "kmernpp.h"
#include "dbgraph.h"
#include "settings.h"
#include "kmer/tstring.h"
#include "dbgraph.h"

using namespace std;

// ============================================================================
// <KMER - NODE-POSITION PAIR> TABLE
// ============================================================================

bool KmerNPPTable::insert(const Kmer& kmer, const NodePosPair& npp)
{
        // choose a representative kmer
        Kmer reprKmer = doubleStranded ? kmer.getRepresentative() : kmer;

        // choose a representative npp
        NodeLength nl = dBG.getSSNode(npp.getNodeID()).getMarginalLength();
        NodePosPair reprNpp = (kmer == reprKmer) ? npp : npp.getRevCompl(nl);

        // insert value in table
        auto it = table.insert(pair<Kmer, NodePosPair>(reprKmer, reprNpp));
        return it.second;
}

NodePosPair KmerNPPTable::find(const Kmer& kmer) const
{
        // choose a representative kmer
        Kmer reprKmer = doubleStranded ? kmer.getRepresentative() : kmer;

        // find the kmer in the table
        auto it = table.find(reprKmer);

        // if it is not found, get out
        if (it == table.end())
                return NodePosPair(0, 0);

        const NodePosPair npp = it->second;
        NodeLength nl = dBG.getSSNode(npp.getNodeID()).getMarginalLength();
        return (kmer == reprKmer) ? npp : npp.getRevCompl(nl);
}
