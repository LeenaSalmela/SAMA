/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


#include <dai/dai_config.h>
#ifdef DAI_WITH_BP


#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <random>
#include <dai/parallelbp.h>
#include <dai/util.h>
#include <dai/properties.h>
#include "dai/exceptions.h"
#include <iomanip>


namespace dai {


using namespace std;


#define DAI_BP_FAST 1
/// \todo Make DAI_BP_FAST a compile-time choice, as it is a memory/speed tradeoff


void ParallelBP::setProperties( const PropertySet &opts ) {
    DAI_ASSERT( opts.hasKey("tol") );
    DAI_ASSERT( opts.hasKey("logdomain") );
    DAI_ASSERT( opts.hasKey("updates") );

    props.tol = opts.getStringAs<Real>("tol");
    props.logdomain = opts.getStringAs<bool>("logdomain");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");

    if( opts.hasKey("maxiter") )
        props.maxiter = opts.getStringAs<size_t>("maxiter");
    else
        props.maxiter = 10000;
    if( opts.hasKey("maxtime") )
        props.maxtime = opts.getStringAs<Real>("maxtime");
    else
        props.maxtime = INFINITY;
    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0;
    if( opts.hasKey("damping") )
        props.damping = opts.getStringAs<Real>("damping");
    else
        props.damping = 0.0;
    if( opts.hasKey("weightdecay") )
            props.weightdecay = opts.getStringAs<bool>("weightdecay");
    else 
            props.weightdecay = false;
    if( opts.hasKey("forcestop") )
            props.forcestop = opts.getStringAs<bool>("forcestop");
    else
            props.forcestop = false;
    if( opts.hasKey("inference") )
        props.inference = opts.getStringAs<Properties::InfType>("inference");
    else
        props.inference = Properties::InfType::SUMPROD;
    if ( opts.hasKey("resinit") )
        props.resinit = opts.getStringAs<Properties::ResInitType>("resinit");
    else 
        props.resinit = Properties::ResInitType::MESSAGE;
    if ( opts.hasKey("modelcount") )
        props.modelcount = opts.getStringAs<size_t>("modelcount");
    else
        props.modelcount = 0;
    if (opts.hasKey("nthreads"))
        props.nthreads = opts.getStringAs<size_t>("nthreads");
    else
        props.nthreads = thread::hardware_concurrency();
    if (opts.hasKey("splashsize"))
        props.splashsize = opts.getStringAs<size_t>("splashsize");
    else
        props.splashsize = 10;
}


PropertySet ParallelBP::getProperties() const {
    PropertySet opts;
    opts.set( "tol", props.tol );
    opts.set( "maxiter", props.maxiter );
    opts.set( "maxtime", props.maxtime );
    opts.set( "verbose", props.verbose );
    opts.set( "logdomain", props.logdomain );
    opts.set( "updates", props.updates );
    opts.set( "damping", props.damping );
    opts.set( "weightdecay", props.weightdecay );
    opts.set( "forcestop", props.forcestop);
    opts.set( "inference", props.inference );
    opts.set( "resinit", props.resinit );
    opts.set( "nthreads", props.nthreads );
    opts.set( "splashsize", props.splashsize);
    return opts;
}


string ParallelBP::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "maxtime=" << props.maxtime << ",";
    s << "verbose=" << props.verbose << ",";
    s << "logdomain=" << props.logdomain << ",";
    s << "updates=" << props.updates << ",";
    if(props.updates == Properties::UpdateType::SPLASH)
            s << "(size=" << props.splashsize << "),"; 
    s << "damping=" << props.damping << ",";
    s << "weightdecay=" << props.weightdecay << ",";
    s << "forcestop=" << props.forcestop << ",";
    s << "inference=" << props.inference << ";";
    s << "resinit=" << props.resinit << ";";
    s << "nthreads=" << props.nthreads << "]";
    return s.str();
}


void ParallelBP::construct() {
    // create edge properties
    _edges.clear();
    _edges.reserve( nrVars() );
    //_eMutex.clear();
    //_eMutex.reserve( nrVars() );
    //_edge2lut.clear();
    _fnode2lut.clear();
    if( props.updates == Properties::UpdateType::SPLASH){
        //_edge2lut.reserve(nrVars());
        _fnode2lut.reserve(nrFactors());
    }
    for( size_t i = 0; i < nrVars(); ++i ) {
        _edges.push_back( vector<EdgeProp>() );
        _edges[i].reserve( nbV(i).size() );
        //_eMutex[i] = std::vector<std::unique_ptr<std::mutex>>(nbV(i).size());
        bforeach( const Neighbor &I, nbV(i) ) {
            EdgeProp newEP;
            newEP.message = Prob( var(i).states() );
            newEP.newMessage = Prob( var(i).states() );

            if( DAI_BP_FAST ) {
                newEP.index.reserve( factor(I).nrStates() );
                for( IndexFor k( var(i), factor(I).vars() ); k.valid(); ++k )
                    newEP.index.push_back( k );
            }

            newEP.residual = 0.0;
            newEP.numUpdate = 0;
            _edges[i].push_back( newEP );
            //std::unique_ptr<std::mutex> mutexptr(new std::mutex());
            //_eMutex[i].emplace_back(mutexptr);
        }
    }
    
    if(props.updates == Properties::UpdateType::SPLASH) {
        for( size_t I = 0; I < nrFactors(); ++I) {
                _fnode2lut.push_back(_vlut.insert(make_pair(std::numeric_limits<double>::max(), I)));
        }    
    } 

    // create old beliefs
    _oldBeliefsV.clear();
    _oldBeliefsV.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i )
        _oldBeliefsV.push_back( Factor( var(i) ) );
    _oldBeliefsF.clear();
    _oldBeliefsF.reserve( nrFactors() );
    for( size_t I = 0; I < nrFactors(); ++I )
        _oldBeliefsF.push_back( Factor( factor(I).vars() ) );
    
    // create update sequence
    _updateSeq.clear();
    _updateSeq.reserve( nrEdges() );
    for( size_t I = 0; I < nrFactors(); I++ )
        bforeach( const Neighbor &i, nbF(I) ){
            _updateSeq.push_back( Edge( i, i.dual ) );
	}
}


void ParallelBP::init() {
    Real c = props.logdomain ? 0.0 : 1.0;
    for( size_t i = 0; i < nrVars(); ++i ) {
        bforeach( const Neighbor &I, nbV(i) ) {
            message( i, I.iter ).fill( c );
            newMessage( i, I.iter ).fill( c );
            if(props.updates == Properties::UpdateType::SPLASH )
                updateResidual( i, I.iter, 0.0 );
        }
    }
    _iters = 0;
}


/*void ParallelBP::findMaxResidual( size_t &i, size_t &_I ) {
    DAI_ASSERT( !_lut.empty() );
    LutType::const_iterator largestEl = _lut.end();
    --largestEl;
    i  = largestEl->second.first;
    _I = largestEl->second.second;
}*/

Real ParallelBP::findMaxNodeResidual(size_t &I, Real& residual) {
        std::lock_guard<std::mutex> guard(vlutMutex);
        DAI_ASSERT( !_vlut.empty() );
        vLutType::const_reverse_iterator largestEl =_vlut.rbegin();
        I = largestEl->second;
        residual = largestEl->first;
        _vlut.erase(_fnode2lut[I]);
        _fnode2lut[I] = _vlut.insert(make_pair(0.0,I));
        return residual;
}


Prob ParallelBP::calcIncomingMessageProduct( size_t I, bool without_i, size_t i ) const {
    Factor Fprod( factor(I) );
    Prob &prod = Fprod.p();
    if( props.logdomain )
        prod.takeLog();

    // Calculate product of incoming messages and factor I
    // TODO what if message is updated by other thread in between?
    bforeach( const Neighbor &j, nbF(I) )
        if( !(without_i && (j == i)) ) {
            // prod_j will be the product of messages coming into j
            Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0 );
            bforeach( const Neighbor &J, nbV(j) )
                if( J != I ) { // for all J in nb(j) \ I
                    std::lock_guard<std::mutex> guard(_edges[j][J.iter].eMutex); // TODO somehow lock all involved edges for the duration of this function?
                    if( props.logdomain )
                        prod_j += message( j, J.iter );
                    else
                        prod_j *= message( j, J.iter );
                }

            // multiply prod with prod_j
            if( !DAI_BP_FAST ) {
                // UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION
                if( props.logdomain )
                    Fprod += Factor( var(j), prod_j );
                else
                    Fprod *= Factor( var(j), prod_j );
            } else {
                // OPTIMIZED VERSION
                size_t _I = j.dual;
                // ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
                const ind_t &ind = index(j, _I);

                for( size_t r = 0; r < prod.size(); ++r )
                    if( props.logdomain )
                        prod.set( r, prod[r] + prod_j[ind[r]] );
                    else
                        prod.set( r, prod[r] * prod_j[ind[r]] );
            }
    }
    return prod;
}


void ParallelBP::calcNewMessage( size_t i, size_t _I ) {
    // calculate updated message I->i
    size_t I = nbV(i,_I);

    Prob marg;
    if( factor(I).vars().size() == 1 ) // optimization
        marg = factor(I).p();
    else {
        Factor Fprod( factor(I) );
        Prob &prod = Fprod.p();
        prod = calcIncomingMessageProduct( I, true, i );

        if( props.logdomain ) {
            prod -= prod.max();
            prod.takeExp();
        }

        // Marginalize onto i
        if( !DAI_BP_FAST ) {
            // UNOPTIMIZED (SIMPLE TO READ, BUT SLOW) VERSION
            if( props.inference == Properties::InfType::SUMPROD )
                marg = Fprod.marginal( var(i) ).p();
            else
                marg = Fprod.maxMarginal( var(i) ).p();
        } else {
            // OPTIMIZED VERSION 
            marg = Prob( var(i).states(), 0.0 );
            // ind is the precalculated IndexFor(i,I) i.e. to x_I == k corresponds x_i == ind[k]
            const ind_t ind = index(i,_I);
            if( props.inference == Properties::InfType::SUMPROD )
                for( size_t r = 0; r < prod.size(); ++r )
                    marg.set( ind[r], marg[ind[r]] + prod[r] );
            else
                for( size_t r = 0; r < prod.size(); ++r )
                    if( prod[r] > marg[ind[r]] )
                        marg.set( ind[r], prod[r] );
            marg.normalize();
        }
    }
    // Store result
    std::lock_guard<std::mutex> guard(_edges[i][_I].eMutex);
    if( props.logdomain )
        newMessage(i,_I) = marg.log();
    else
        newMessage(i,_I) = marg;
    
    // Update the residual if necessary
    if(props.updates == Properties::UpdateType::SPLASH){
            if( props.weightdecay)
                    updateResidual( i, _I , dist( newMessage( i, _I ), message( i, _I ), DISTLINF ) / (Real) std::max((size_t) 1, numUpdates(i,_I) ));
            else
                    updateResidual( i, _I , dist( newMessage( i, _I ), message( i, _I ), DISTLINF ));
    }
}


// BP::run does not check for NANs for performance reasons
// Somehow NaNs do not often occur in BP...
Real ParallelBP::run() {
    if( props.verbose >= 1 )
        cerr << "Starting " << identify() << "...";
    //if( props.verbose >= 3)
        cerr << endl;

    double tic = toc();
    
    std::ofstream outV;
    std::ofstream nonConvVar;
    std::ofstream outF;
    std::ofstream outE;
    
    size_t numVarNonConv_old = nrVars();
    size_t timesNVNCincreased = 0;
    
    if( props.verbose >=3){
        string filenameV = name()+"_"+to_string((int)props.inference)+"_"+to_string((int)props.updates)+"_"+to_string((int)props.modelcount)+"_variableBeliefDiff.tsv";
        string filenameF = name()+"_"+to_string((int)props.inference)+"_"+to_string((int)props.updates)+"_"+to_string((int)props.modelcount)+"_factorBeliefDiff.tsv";
        string filenameNC = name()+"_"+to_string((int)props.inference)+"_"+to_string((int)props.updates)+"_"+to_string((int)props.modelcount)+"_notConverged.tsv";
        string filenameE = name()+"_"+to_string((int)props.inference)+"_"+to_string((int)props.updates)+"_"+to_string((int)props.modelcount)+"_messageFreq.tsv";
        outV.open(filenameV);
        outF.open(filenameF);
        nonConvVar.open(filenameNC);
        outE.open(filenameE);
    }

    // do several passes over the network until maximum number of iterations has
    // been reached or until the maximum belief difference is smaller than tolerance
    Real maxDiff = INFINITY;
    for( ; _iters < props.maxiter && maxDiff > props.tol && (toc() - tic) < props.maxtime; _iters++ ) {
        if ( props.updates == Properties::UpdateType::SPLASH) {
                for( size_t i = 0; i < nrVars(); ++i )
                        bforeach( const Neighbor &I, nbV(i) ){
                                calcNewMessage( i, I.iter );
                                // Just sending them all once more to see if maxdiff also shows convergence (might remove this for efficiency gain)
                                if (_iters == 1) 
                                        updateMessage(i, I.iter);
                        }
                //initialise vertex residuals (try to make sure splashes will be evenly spaced)
                if (_iters == 0){
                        for(size_t I = 0; I < nrFactors(); I++)
                                getCurrNodeResidual(I);
                        Real maxRes = getMaxNodeResidual();
                        std::vector<size_t> Fids(nrFactors());
                        std::iota(Fids.begin(), Fids.end(), 1);
                        std::shuffle(Fids.begin(), Fids.end(), std::mt19937{std::random_device{}()});
                        for(size_t I = 0; I < nrFactors(); I++){
                                updateNodeResidual(I, maxRes + Fids[I]);
                        }
                }
            std::vector<size_t> splashCounters(props.nthreads, 0);
            std::vector<std::thread> wt(props.nthreads);
            for (size_t i = 0; i < wt.size(); i++)
                    wt[i] = std::thread(&ParallelBP::splashThreadWork, this, i, tic,  ref(splashCounters[i]));
            
            for_each(wt.begin(), wt.end(), mem_fn(&thread::join));
            
            cerr << name() << ":: sent " << std::accumulate(splashCounters.begin(), splashCounters.end(), 0) << " splashes" << endl;
                
        } else if( props.updates == Properties::UpdateType::PARALL ) {
            // Parallel updates
            for( size_t i = 0; i < nrVars(); ++i )
                bforeach( const Neighbor &I, nbV(i) )
                    calcNewMessage( i, I.iter );

            for( size_t i = 0; i < nrVars(); ++i )
                bforeach( const Neighbor &I, nbV(i) )
                    updateMessage( i, I.iter );
        } else {
            // Sequential updates
            random_shuffle( _updateSeq.begin(), _updateSeq.end(), rnd );

            bforeach( const Edge &e, _updateSeq ) {
                calcNewMessage( e.first, e.second );
                updateMessage( e.first, e.second );
            }
        }

        // calculate new beliefs and compare with old ones
        maxDiff = -INFINITY;
        int numNotC = 0;
        for( size_t i = 0; i < nrVars(); ++i ) {
            Factor b( beliefV(i) );
            if (props.verbose >= 3) {
                outV << dist( b, _oldBeliefsV[i], DISTLINF ) << "\t";
                if (dist( b, _oldBeliefsV[i], DISTLINF ) > props.tol) {
                    nonConvVar << b.vars() << "\t";
                    numNotC++;
                }
            }
            maxDiff = std::max( maxDiff, dist( b, _oldBeliefsV[i], DISTLINF ) );
            _oldBeliefsV[i] = b;
        }
        if (props.verbose >= 3) {
            outV << endl;
            nonConvVar << endl;
        }
        for( size_t I = 0; I < nrFactors(); ++I ) {
            Factor b( beliefF(I) );
            if (props.verbose >= 3)
                outF << dist( b, _oldBeliefsF[I], DISTLINF ) << "\t";
            maxDiff = std::max( maxDiff, dist( b, _oldBeliefsF[I], DISTLINF ) );
            _oldBeliefsF[I] = b;
        }
        
        if (props.verbose >= 3)
            outF << endl;

        if( props.verbose >= 3 ) {
                cerr << name() << "::run: maxdiff " << maxDiff << ", after " << _iters+1 << " passes, " << numNotC << " variable(s) not converged (" << toc() - tic << " seconds)." << endl;}
        
        if( props.forcestop ){
                if (numNotC >= numVarNonConv_old)
                        timesNVNCincreased += 1;
                if (timesNVNCincreased == 3) // TODO make this less arbitrary / a parameter
                        break;
                numVarNonConv_old = numNotC;
        }
    }
    cout.precision(2);
    if( maxDiff > _maxdiff )
        _maxdiff = maxDiff;
    
    if (props.verbose >= 3){
            for (vector<EdgeProp> EV : _edges){
                    for(EdgeProp ep: EV)
                            outE << ep.numUpdate << "\n";
            }
    }

    if( props.verbose >= 1 ) {
        if( maxDiff > props.tol ) {
            if( props.verbose == 1 )
                cerr << endl;
                cerr << name() << "::run:  WARNING: not converged after " << _iters << " passes (" << toc() - tic << " seconds)...final maxdiff:" << maxDiff << endl;
        } else {
            if( props.verbose >= 3 )
                cerr << name() << "::run:  ";
                cerr << "converged in " << _iters << " passes (" << toc() - tic << " seconds)." << endl;
        }
    }

    return maxDiff;
}


void ParallelBP::calcBeliefV( size_t i, Prob &p ) const { // TODO not used during a parallel step?
    p = Prob( var(i).states(), props.logdomain ? 0.0 : 1.0 );
    bforeach( const Neighbor &I, nbV(i) )
        if( props.logdomain )
            p += newMessage( i, I.iter );
        else
            p *= newMessage( i, I.iter );
}


Factor ParallelBP::beliefV( size_t i ) const {
    Prob p;
    calcBeliefV( i, p );

    if( props.logdomain ) {
        p -= p.max();
        p.takeExp();
    }
    p.normalize();

    return( Factor( var(i), p ) );
}


Factor ParallelBP::beliefF( size_t I ) const {
    Prob p;
    calcBeliefF( I, p );

    if( props.logdomain ) {
        p -= p.max();
        p.takeExp();
    }
    p.normalize();

    return( Factor( factor(I).vars(), p ) );
}


vector<Factor> ParallelBP::beliefs() const {
    vector<Factor> result;
    for( size_t i = 0; i < nrVars(); ++i )
        result.push_back( beliefV(i) );
    for( size_t I = 0; I < nrFactors(); ++I )
        result.push_back( beliefF(I) );
    return result;
}


Factor ParallelBP::belief( const VarSet &ns ) const {
    if( ns.size() == 0 )
        return Factor();
    else if( ns.size() == 1 )
        return beliefV( findVar( *(ns.begin() ) ) );
    else {
        size_t I;
        for( I = 0; I < nrFactors(); I++ )
            if( factor(I).vars() >> ns )
                break;
        if( I == nrFactors() )
            DAI_THROW(BELIEF_NOT_AVAILABLE);
        return beliefF(I).marginal(ns);
    }
}


Real ParallelBP::logZ() const {
    Real sum = 0.0;
    for( size_t i = 0; i < nrVars(); ++i )
        sum += (1.0 - nbV(i).size()) * beliefV(i).entropy();
    for( size_t I = 0; I < nrFactors(); ++I )
        sum -= dist( beliefF(I), factor(I), DISTKL );
    return sum;
}


void ParallelBP::init( const VarSet &ns ) {
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); ++n ) {
        size_t ni = findVar( *n );
        bforeach( const Neighbor &I, nbV( ni ) ) {
            Real val = props.logdomain ? 0.0 : 1.0;
            message( ni, I.iter ).fill( val );
            newMessage( ni, I.iter ).fill( val );
            if( props.updates == Properties::UpdateType::SPLASH )
                updateResidual( ni, I.iter, 0.0 );
        }
    }
    _iters = 0;
}

void ParallelBP::getCurrNodeResidual(size_t I) {
        Real res = 0.0;
        bforeach(const Neighbor &j, nbF(I)){
                size_t _I = j.dual;
                bforeach( const Neighbor &J, nbV(j) ){
                        size_t _J = J.iter;
                        if( _J!=_I ){
                                std::lock_guard<std::mutex> guard(_edges[j][_J].eMutex);
                                res = std::max(res, residual(j, _J));
                        }
                }
        }
        updateNodeResidual(I, res);
}

void ParallelBP::calcCurrNodeResidual(size_t I) {
        Real new_res = 0.0;
        bforeach(const Neighbor &j, nbF(I)){
                Prob prod_j_last( var(j).states(), props.logdomain ? 0.0 : 1.0 );
                Prob prod_j_next( var(j).states(), props.logdomain ? 0.0 : 1.0 );
                bforeach( const Neighbor &J, nbV(j) )
                if( J != I ) {
                        std::lock_guard<std::mutex> guard(_edges[j][J.iter].eMutex); // for all J in nb(j) \ I
                        if( props.logdomain ) {
                                prod_j_last += message( j, J.iter );
                                prod_j_next += newMessage( j, J.iter );
                        } else {
                                prod_j_last*= message( j, J.iter );
                                prod_j_next*= newMessage( j, J.iter );
                        }
                }
                new_res = std::max(dist(prod_j_next, prod_j_last, DISTLINF), new_res);
        }
        updateNodeResidual(I, new_res);
}

void ParallelBP::updateMessage( size_t i, size_t _I ) {
    std::lock_guard<std::mutex> guard(_edges[i][_I].eMutex);
    if( recordSentMessages )
        _sentMessages.push_back(make_pair(i,_I));
    if( props.damping == 0.0 ) {
        message(i,_I) = newMessage(i,_I);
        if(props.updates == Properties::UpdateType::SPLASH)
            updateResidual( i, _I, 0.0 );
    } else {
        if( props.logdomain )
            message(i,_I) = (message(i,_I) * props.damping) + (newMessage(i,_I) * (1.0 - props.damping));
        else
            message(i,_I) = (message(i,_I) ^ props.damping) * (newMessage(i,_I) ^ (1.0 - props.damping));
        if( props.updates == Properties::UpdateType::SPLASH)
            updateResidual( i, _I, dist( newMessage(i,_I), message(i,_I), DISTLINF ) );
    }
    addUpdate(i,_I);
}


void ParallelBP::updateResidual( size_t i, size_t _I, Real r ) { //TODO always locked edge before use?
        EdgeProp* pEdge = &_edges[i][_I];
        pEdge->residual = r;

    // rearrange look-up table (delete and reinsert new key)
    //_lut.erase( _edge2lut[i][_I] );
    //_edge2lut[i][_I] = _lut.insert( make_pair( r, make_pair(i, _I) ) );
}

void ParallelBP::updateNodeResidual(size_t I, Real r) {
        std::lock_guard<std::mutex> guard(vlutMutex);
        _vlut.erase(_fnode2lut[I]);
        _fnode2lut[I] = _vlut.insert(make_pair(r,I));
}

/*void ParallelBP::increasePriority(size_t i, size_t _I, Real r) {
    EdgeProp* pEdge = &_edges[i][_I];
    pEdge->residual += r;

    // rearrange look-up table (delete and reinsert new key)
    _lut.erase( _edge2lut[i][_I] );
    _edge2lut[i][_I] = _lut.insert( make_pair( pEdge->residual, make_pair(i, _I) ) );
}*/


void ParallelBP::splashThreadWork(const size_t tID, const double tic, size_t& splashCounter) {
        // Get highest residual root
        size_t I_root; Real I_residual;
        while(findMaxNodeResidual(I_root, I_residual) > props.tol && (toc() - tic) < props.maxtime){
                // Get splash
                // TODO fix/parameterise search depth
                std::vector<size_t> sendOrder = getSplash(I_root, props.splashsize);
                //std::unique_lock<std::mutex> vlut_guard(vlutMutex);
                //std::vector<size_t> sendOrder = getResidualSplash(I_root, 30, _fnode2lut, props.tol); //TODO avoid root_resid = 0 bcs ResidualSplash doesn't work anymore
                //vlut_guard.unlock();
                std::set<size_t> influencedNodes(sendOrder.begin(), sendOrder.end());
                // Pass messages up and down splash
                if (props.verbose > 2){
                std::lock_guard<std::mutex> lg(writeMutex);
                cerr << "T" << tID << std::setprecision(5) << ": Splash root " << I_root << ", root node residual: " << I_residual << ", Splash size: " << sendOrder.size() << "\n";
                }
                splashCounter += 1;
                //cout << "Downward pass" << endl;
                for (size_t idx = sendOrder.size(); --idx > 0 ;){
                        size_t I = sendOrder[idx];
                        bforeach(const Neighbor &i, nbF(I)){
                                size_t _I = i.dual;
                                //calcNewMessage(i, _I);
                                updateMessage(i, _I);
                                //cout << I << "\t" << i.node << " (\t";
                                // I->i has been updated, which means that residuals for all
                                // J->j with J in nb[i]\I and j in nb[J]\i have to be updated
                                bforeach( const Neighbor &J, nbV(i) ) {
                                        if( J.iter != _I ) {
                                                bforeach( const Neighbor &j, nbF(J) ) {
                                                        size_t _J = j.dual;
                                                        if( j != i ) {
                                                                calcNewMessage( j, _J );
                                                        }
                                                }
                                        }
                                }
                                //cout << ')' << endl;
                        }
                }
                //cout << "Upward pass" << endl;
                for (size_t idx = 0; idx < sendOrder.size(); idx++){
                        size_t I = sendOrder[idx];
                        bforeach(const Neighbor &i, nbF(I)){
                                size_t _I = i.dual;
                                //calcNewMessage(i, _I);
                                updateMessage(i, _I);
                                //cout << I << "\t" << i.node << " (\t";
                                // I->i has been updated, which means that residuals for all
                                // J->j with J in nb[i]\I and j in nb[J]\i have to be updated
                                bforeach( const Neighbor &J, nbV(i) ) {
                                        if( J.iter != _I ) {
                                                bforeach( const Neighbor &j, nbF(J) ) {
                                                        size_t _J = j.dual;
                                                        if( j != i ) {
                                                                calcNewMessage( j, _J );
                                                                influencedNodes.insert(J.node);
                                                                bforeach(const Neighbor &K, nbV(j))
                                                                if( K.node != J.node)
                                                                        influencedNodes.insert(K);
                                                                //cout << J.node << "  " << j.node << "\t";
                                                        }
                                                }
                                        }
                                }
                                //cout << ')' << endl;
                        }
                }
                // Update Node residuals (propagate influence outside splash)
                for (size_t I: influencedNodes){
                        //for(size_t I = 0; I < nrFactors(); I++){
                        getCurrNodeResidual(I);
                        //cout << _fnode2lut[I]->first << "\t";
                }
                //cout << endl;
        } 
}


} // end of namespace dai


#endif
