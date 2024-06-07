/*  This file is part of libDAI - http://www.libdai.org/
 *
 *  Copyright (c) 2006-2011, The libDAI authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


/// \file
/// \brief Defines class BP, which implements (Loopy) Belief Propagation
/// \todo Consider using a priority_queue for maximum residual schedule


#ifndef __defined_libdai_parallelbp_h
#define __defined_libdai_parallelbp_h


#include <dai/dai_config.h>
#ifdef DAI_WITH_PARALLELBP


#include <string>
#include <mutex>
#include <thread>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/properties.h>
#include <dai/enum.h>
#include <dai/bp.h>


namespace dai {


/// Approximate inference algorithm "(Loopy) Belief Propagation"
/** The Loopy Belief Propagation algorithm uses message passing
 *  to approximate marginal probability distributions ("beliefs") for variables
 *  and factors (more precisely, for the subset of variables depending on the factor).
 *  There are two variants, the sum-product algorithm (corresponding to 
 *  finite temperature) and the max-product algorithm (corresponding to 
 *  zero temperature).
 *
 *  The messages \f$m_{I\to i}(x_i)\f$ are passed from factors \f$I\f$ to variables \f$i\f$. 
 *  In case of the sum-product algorith, the update equation is: 
 *    \f[ m_{I\to i}(x_i) \propto \sum_{x_{N_I\setminus\{i\}}} f_I(x_I) \prod_{j\in N_I\setminus\{i\}} \prod_{J\in N_j\setminus\{I\}} m_{J\to j}\f]
 *  and in case of the max-product algorithm:
 *    \f[ m_{I\to i}(x_i) \propto \max_{x_{N_I\setminus\{i\}}} f_I(x_I) \prod_{j\in N_I\setminus\{i\}} \prod_{J\in N_j\setminus\{I\}} m_{J\to j}\f]
 *  In order to improve convergence, the updates can be damped. For improved numerical stability,
 *  the updates can be done in the log-domain alternatively.
 *
 *  After convergence, the variable beliefs are calculated by:
 *    \f[ b_i(x_i) \propto \prod_{I\in N_i} m_{I\to i}(x_i)\f]
 *  and the factor beliefs are calculated by:
 *    \f[ b_I(x_I) \propto f_I(x_I) \prod_{j\in N_I} \prod_{J\in N_j\setminus\{I\}} m_{J\to j}(x_j) \f]
 *  The logarithm of the partition sum is calculated by:
 *    \f[ \log Z = \sum_i (1 - |N_i|) \sum_{x_i} b_i(x_i) \log b_i(x_i) - \sum_I \sum_{x_I} b_I(x_I) \log \frac{b_I(x_I)}{f_I(x_I)} \f]
 *
 *  For the max-product algorithm, a heuristic way of finding the MAP state (the 
 *  joint configuration of all variables which has maximum probability) is provided
 *  by the findMaximum() method, which can be called after convergence.
 *
 *  \note There are two implementations, an optimized one (the default) which caches IndexFor objects,
 *  and a slower, less complicated one which is easier to maintain/understand. The slower one can be 
 *  enabled by defining DAI_BP_FAST as false in the source file.
 */
class ParallelBP : public DAIAlgFG {
    protected:
        /// Type used for index cache
        typedef std::vector<size_t> ind_t;
        /// Type used for storing edge properties
        struct EdgeProp {
            EdgeProp(): message(), newMessage(), residual(), numUpdate() {}
            EdgeProp(const EdgeProp& ep) {
                    std::lock_guard<std::mutex> lg(ep.eMutex);
                    index = ep.index;
                    message = ep.message;
                    newMessage = ep.newMessage;
                    residual = ep.residual;
                    numUpdate = ep.numUpdate;
            }
            EdgeProp(EdgeProp&& ep){ // std::exchange since C++14
                    std::lock_guard<std::mutex> lg(ep.eMutex);
                    std::swap(index, ep.index);
                    std::swap(message, ep.message);
                    std::swap(newMessage, ep.newMessage);
                    std::swap(residual, ep.residual);
                    std::swap(numUpdate, ep.numUpdate);
            }
            EdgeProp& operator=(EdgeProp&& ep){
                    if(&ep != this){
                        std::lock_guard<std::mutex> lg1(this->eMutex), lg2(ep.eMutex);
                        std::swap(index, ep.index);
                        std::swap(message, ep.message);
                        std::swap(newMessage, ep.newMessage);
                        std::swap(residual, ep.residual);
                        std::swap(numUpdate, ep.numUpdate);
                    }
                    return *this;
            }
            EdgeProp& operator=(const EdgeProp& ep){
                    if (&ep != this){
                        std::lock_guard<std::mutex> lg1(this->eMutex), lg2(ep.eMutex);
                        index = ep.index;
                        message = ep.message;
                        newMessage = ep.newMessage;
                        residual = ep.residual;
                        numUpdate = ep.numUpdate;
                    }
                    return *this;
            }
            /// Index cached for this edge
            ind_t  index;
            /// Old message living on this edge
            Prob   message;
            /// New message living on this edge
            Prob   newMessage;
            /// Residual for this edge
            Real   residual;
            /// Count number of updates (For Weight decay BP)
            size_t numUpdate;
            /// Mutex for this edge
            mutable std::mutex eMutex; //avoid default move/copy/assignment of EdgeProp, bcs deleted in mutex
        };
        /// Stores all edge properties
        std::vector<std::vector<EdgeProp> > _edges;
        /// One mutex for each edge
        //mutable std::vector<std::vector<std::unique_ptr<std::mutex>>> _eMutex;
        /// Type of lookup table (only used for maximum-residual BP) (priority queue, but with easier insertion and removal)
        //typedef std::multimap<Real, std::pair<size_t, size_t> > LutType;
        /// Lookup table (only used for maximum-residual BP)
        // std::vector<std::vector<LutType::iterator> > _edge2lut;
        /// Lookup table (only used for maximum-residual BP)
        // LutType _lut;
        /// vertex residual priority queue (we only consider factor vertices here)
        typedef std::multimap<Real,size_t> vLutType;
        vLutType _vlut;
        std::vector<vLutType::iterator> _fnode2lut;
        /// Maximum difference between variable beliefs encountered so far
        Real _maxdiff;
        /// Number of iterations needed
        size_t _iters;
        /// The history of message updates (only recorded if \a recordSentMessages is \c true)
        std::vector<std::pair<size_t, size_t> > _sentMessages;
        /// Stores variable beliefs of previous iteration
        std::vector<Factor> _oldBeliefsV;
        /// Stores factor beliefs of previous iteration
        std::vector<Factor> _oldBeliefsF;
        /// Stores the update schedule
        std::vector<Edge> _updateSeq;
        
        //mutable std::mutex edgesMutex;
        mutable std::mutex writeMutex;
        mutable std::mutex vlutMutex;
        
        //std::atomic<size_t> splashCounter;

    public:
        /// Parameters for BP
        struct Properties {
            /// Enumeration of possible update schedules
            /** The following update schedules have been defined:
             *  - PARALL parallel updates
             *  - SEQFIX sequential updates using a fixed sequence
             *  - SEQRND sequential updates using a random sequence
             *  - SEQMAX maximum-residual updates [\ref EMK06]
             */
            DAI_ENUM(UpdateType,PARALL,SPLASH);

            /// Enumeration of inference variants
            /** There are two inference variants:
             *  - SUMPROD Sum-Product
             *  - MAXPROD Max-Product (equivalent to Min-Sum)
             */
            DAI_ENUM(InfType,SUMPROD,MAXPROD);
            
            /// Type of residual initialisation (only used for SEQMAX0L)
            DAI_ENUM(ResInitType,MESSAGE,UNIFORM);

            /// Verbosity (amount of output sent to stderr)
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Maximum time (in seconds)
            double maxtime;

            /// Tolerance for convergence test
            Real tol;

            /// Whether updates should be done in logarithmic domain or not
            bool logdomain;
            
            /// Perform weight decay residual updates (Only for SEQMAX)?
            bool weightdecay;
            
            /// Force a stop of LBP based on how many variables have converged
            bool forcestop;

            /// Damping constant (0.0 means no damping, 1.0 is maximum damping)
            Real damping;
            
            /// Number of threads to use
            size_t nthreads;
            
            /// Splash size
            size_t splashsize;

            /// Message update schedule
            UpdateType updates;

            /// Inference variant
            InfType inference;
            
            /// Initialisation variant
            ResInitType resinit;
            
            /// For output naming
            size_t modelcount;
        } props;

        /// Specifies whether the history of message updates should be recorded
        bool recordSentMessages;

    public:
    /// \name Constructors/destructors
    //@{
        /// Default constructor
        ParallelBP() : DAIAlgFG(), _edges(), _vlut(), _maxdiff(0.0), _iters(0U), _sentMessages(), _oldBeliefsV(), _oldBeliefsF(), _updateSeq(), props(), recordSentMessages(false) {}

        /// Construct from FactorGraph \a fg and PropertySet \a opts
        /** \param fg Factor graph.
         *  \param opts Parameters @see Properties
         */
        ParallelBP( const FactorGraph & fg, const PropertySet &opts ) : DAIAlgFG(fg), _edges(), _maxdiff(0.0), _iters(0U), _sentMessages(), _oldBeliefsV(), _oldBeliefsF(), _updateSeq(), props(), recordSentMessages(false) {
            setProperties( opts );
            construct();
        }

        /// Copy constructor
        ParallelBP( const ParallelBP &x ) : DAIAlgFG(x), _edges(x._edges), _vlut(x._vlut), _fnode2lut(x._fnode2lut), _maxdiff(x._maxdiff), _iters(x._iters), _sentMessages(x._sentMessages), _oldBeliefsV(x._oldBeliefsV), _oldBeliefsF(x._oldBeliefsF), _updateSeq(x._updateSeq), props(x.props), recordSentMessages(x.recordSentMessages) {
            for(vLutType::iterator l = _vlut.begin(); l != _vlut.end(); ++l )
                    _fnode2lut[l->second] = l;
        }

        /// Assignment operator
        ParallelBP& operator=( const ParallelBP &x ) {
            if( this != &x ) {
                DAIAlgFG::operator=( x );
                _edges = x._edges;
                _vlut = x._vlut;
                _fnode2lut = x._fnode2lut;
                for(vLutType::iterator l = _vlut.begin(); l != _vlut.end(); ++l )
                        _fnode2lut[l->second] = l;
                _maxdiff = x._maxdiff;
                _iters = x._iters;
                _sentMessages = x._sentMessages;
                _oldBeliefsV = x._oldBeliefsV;
                _oldBeliefsF = x._oldBeliefsF;
                _updateSeq = x._updateSeq;
                props = x.props;
                recordSentMessages = x.recordSentMessages;
            }
            return *this;
        }
    //@}

    /// \name General InfAlg interface
    //@{
        virtual ParallelBP* clone() const { return new ParallelBP(*this); }
        virtual ParallelBP* construct( const FactorGraph &fg, const PropertySet &opts ) const { return new ParallelBP( fg, opts ); }
        virtual std::string name() const { return "PARALLELBP"; }
        virtual Factor belief( const Var &v ) const { return beliefV( findVar( v ) ); }
        virtual Factor belief( const VarSet &vs ) const;
        virtual Factor beliefV( size_t i ) const;
        virtual Factor beliefF( size_t I ) const;
        virtual std::vector<Factor> beliefs() const;
        virtual Real logZ() const;
        /** \pre Assumes that run() has been called and that \a props.inference == \c MAXPROD
         */
        std::vector<size_t> findMaximum() const { return dai::findMaximum( *this ); }
        virtual void init();
        virtual void init( const VarSet &ns );
        virtual Real run();
        virtual Real maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
        virtual void setMaxIter( size_t maxiter ) { props.maxiter = maxiter; }
        virtual void setProperties( const PropertySet &opts );
        virtual PropertySet getProperties() const;
        virtual std::string printProperties() const;
    //@}

    /// \name Additional interface specific for BP
    //@{
        /// Returns history of which messages have been updated
        const std::vector<std::pair<size_t, size_t> >& getSentMessages() const {
            return _sentMessages;
        }

        /// Clears history of which messages have been updated
        void clearSentMessages() { _sentMessages.clear(); }
    //@}

    protected:
        /// Returns constant reference to message from the \a _I 'th neighbor of variable \a i to variable \a i
        const Prob & message(size_t i, size_t _I) const { return _edges[i][_I].message; }
        /// Returns reference to message from the \a _I 'th neighbor of variable \a i to variable \a i
        Prob & message(size_t i, size_t _I) { return _edges[i][_I].message; }
        void setMessage(size_t i, size_t _I , Prob m) {std::lock_guard<std::mutex> guard(_edges[i][_I].eMutex); _edges[i][_I].message = m;}
        /// Returns constant reference to updated message from the \a _I 'th neighbor of variable \a i to variable \a i
        const Prob & newMessage(size_t i, size_t _I) const { return _edges[i][_I].newMessage; }
        /// Returns reference to updated message from the \a _I 'th neighbor of variable \a i to variable \a i
        Prob & newMessage(size_t i, size_t _I) { return _edges[i][_I].newMessage; }
        void setNewMessage(size_t i, size_t _I , Prob m) {std::lock_guard<std::mutex> guard(_edges[i][_I].eMutex); _edges[i][_I].newMessage = m;}
        /// Returns constant reference to cached index for the edge between variable \a i and its \a _I 'th neighbor
        const ind_t & index(size_t i, size_t _I) const { return _edges[i][_I].index; }
        /// Returns reference to cached index for the edge between variable \a i and its \a _I 'th neighbor
        ind_t & index(size_t i, size_t _I) { return _edges[i][_I].index; }
        /// Returns constant reference to residual for the edge between variable \a i and its \a _I 'th neighbor
        const Real & residual(size_t i, size_t _I) const { return _edges[i][_I].residual; }
        /// Returns reference to residual for the edge between variable \a i and its \a _I 'th neighbor
        Real & residual(size_t i, size_t _I) { return _edges[i][_I].residual; }
        /// Update number of updates
        void addUpdate(size_t i, size_t _I) { _edges[i][_I].numUpdate += 1; }
        const size_t numUpdates(size_t i, size_t _I) const { return _edges[i][_I].numUpdate; }

        /// Calculate the product of factor \a I and the incoming messages
        /** If \a without_i == \c true, the message coming from variable \a i is omitted from the product
         *  \note This function is used by calcNewMessage() and calcBeliefF()
         */
        virtual Prob calcIncomingMessageProduct( size_t I, bool without_i, size_t i ) const;
        /// Calculate the updated message from the \a _I 'th neighbor of variable \a i to variable \a i
        virtual void calcNewMessage( size_t i, size_t _I );
        /// Replace the "old" message from the \a _I 'th neighbor of variable \a i to variable \a i by the "new" (updated) message
        void updateMessage( size_t i, size_t _I );
        /// Get the current value of the node residual based on all messages
        void getCurrNodeResidual(size_t I);
        void calcCurrNodeResidual(size_t I);
        /// Set the residual (difference between new and old message) for the edge between variable \a i and its \a _I 'th neighbor to \a r
        void updateResidual( size_t i, size_t _I, Real r );
        /// Set the node residual (i.e. sup_nb(I) r_m(I->nb(I))) for factor node I to r
        void updateNodeResidual(size_t I, Real r);
	/// increase the priority for SEQMAX0L (saved as residual object!!) for the edge between variable \a i and its \a _I 'th neighbor to \a r
	// void increasePriority(size_t i, size_t _I, Real r);
        /// Finds the edge which has the maximum residual (difference between new and old message)
        // void findMaxResidual( size_t &i, size_t &_I );
        /// Finds the Factor node that has the highest node residual (See Gonzalez et al. 2009)
        Real findMaxNodeResidual(size_t &I, Real& residual);
        Real getMaxNodeResidual(){
                std::lock_guard<std::mutex> guard(vlutMutex);
                DAI_ASSERT( !_vlut.empty() );
                vLutType::const_iterator largestEl =_vlut.end();
                -- largestEl;
                return largestEl->first;
        }
        /// Calculates unnormalized belief of variable \a i
        virtual void calcBeliefV( size_t i, Prob &p ) const;
        /// Calculates unnormalized belief of factor \a I
        virtual void calcBeliefF( size_t I, Prob &p ) const {
            p = calcIncomingMessageProduct( I, false, 0 );
        }

        /// Helper function for constructors
        virtual void construct();
        
        void splashThreadWork(const size_t tID, const double tic, size_t& splashCounter);
};


} // end of namespace dai


#endif


#endif
