#ifndef DD_H
#define DD_H

#include "typedef.h"
#include <mpi.h>

class Subdomain;
class RHS;

//! The abstract class for a domain decomposition (the core part of the preconditioner)
/*!
 * This class is overloaded by specific domain decomposition routines.
 */
class DD
{
  public:
    //! The MPI communicator for the DD
    MPI_Comm mComm;

  public:
    //! The standard constructor
    DD() {};

    //! The copy constructor
    DD( const DD &other ): mComm(other.mComm) {};

    //! The destructor
    virtual ~DD() {};

    //! The copy function to get a pointer to a copy
    virtual DD* copy() const = 0;

    //! A function that returns the communicator
    /*!
     * \return the MPI communicator
     */
    MPI_Comm getComm() const { return mComm; }

    //! A function that factorizes all local linear systems
    virtual void factorize() = 0;

    //! A function that assembles the RHS
    /*!
     * \param rhs The RHS used for assembling 
     * \param res The result given as a vector of complex vectors. Each entry in the
     *        vector represents one subdomain
     */
    virtual void assembleRHS( const RHS *rhs, vector<ComplexVector> &res ) const = 0;

    //! A function that applies the preconditioner 
    /*!
     * \param x The vector that the preconditioner should be applied to
     *          This vector is assumed to be given as a vector of complex vectors for
     *          which each entry corresponds to the local vector in each subdomain.
     * \param res The result given as a vector of complex vectors. Each entry in the
     *        vector corresponds to one subdomain
     */
    virtual void applyPrecond( const vector<ComplexVector> &x, vector<ComplexVector> &res ) const = 0;

    //! A function that applies the system matrix
    /*!
     * \param x The vector that the system matrix should be applie to.
     *          This vector is assumed to be given as a vector of complex vectors for
     *          which each entry corresponds to the local vector in each subdomain.
     * \param res The result given as a vector of complex vectors. Each entry in the
     *        vector corresponds to one subdomain
     */
    virtual void applySystemMatrix( const vector<ComplexVector> &x, vector<ComplexVector> &res ) const = 0;

    //! A function that returns the dimension of the underlying problem
    virtual UInt getDim() const = 0;

    //! A function that returns the global mesh size
    virtual Double getMeshSize() const = 0;
    
    //! A function that expands the values in each subdomain
    /*!
     * The input is assumed to be given as a vector of vector of values corresponding to all  interior DOFs
     * These DOFs are called  \fs\boldsymbol{\Omega}_{ij}\fs in the paper
     *
     * \p x The input values
     * \p res The output values, is overwritten
     * \p shiftedLR The flag if the values correspond to ones shifted in the Left-Right direction
     * \p shiftedUD The flag if the values correspond to ones shifted in the Up-Down direction
     */
    virtual void expand( const vector<ComplexVector> &x, vector<ComplexVector> &res, bool shiftedLR = false, bool shiftedUD=false ) const = 0;
   
    //! A function that communicates trace information from top right to bottom left
    /*!
     * \param xEx The input and output vector of values in each subdomain, the values are assumed to be given in expanded form \sa expand.
     * \param shifted The flag if the input/output vector holds the values in shifted form (or not)
     *
     * This function takes the bottom/left trace (traceN) and communicates them to the bottom/left element.
     */
    virtual void communicateBL2TRTraces( vector<ComplexVector> &xEx, bool shifted=false ) const = 0;
    
    //! The function that communicates trace information from top right to bottom left
    /*!
     * \param xEx The input and output vector of values in each subdomain, the values are assumed to be given in expanded form \sa expand.
     * \param shifted The flag if the input/output vector holds the values in shifted form (or not)
     *
     * This function takes the bottom/left trace (traceN) and communicates them to the bottom/left element.
     */
    virtual void communicateTR2BLTraces( vector<ComplexVector> &xEx, bool shifted=false ) const = 0;

    //! A funciton that returns the total number of cells
    virtual UInt getNCells() const = 0;

    //! A function that returns a pointer to a cell
    /*!
     *  \p i The cell index.
     */
    virtual Subdomain* getCell( UInt i ) const = 0;
};

#endif
