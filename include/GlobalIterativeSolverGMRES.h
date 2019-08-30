#ifndef GLOBALITERATIVESOLVERGMRES_H
#define GLOBALITERATIVESOLVERGMRES_H

#include "GlobalIterativeSolver.h"
#include <mpi.h>

//! A class that implements a (preconditioned) GMRES method for the global solution of a problem defined in a domain decomposition.
/*!
 * This class uses all the functions provided by the domain decomposition and solves
 * the global system in parallel using a preconditioned GMRES method.
 */
class GlobalIterativeSolverGMRES : public GlobalIterativeSolver
{
  protected:
    //! The tolerance for the GMRES method
    Double mTol;

    //! The maximum number of iterations used in the GMRES method
    UInt mMaxIter;

    //! A flag if the GMRES method should print out residual info in each iteration
    bool mVerbose;

    //! The vector that stores all residual information
    mutable vector<Double> mLog;

   
  public:
    //! The standard constructor
    GlobalIterativeSolverGMRES();

    //! The copy constructor
    GlobalIterativeSolverGMRES( const GlobalIterativeSolverGMRES &other );

    //! The destructor
    virtual ~GlobalIterativeSolverGMRES() {};

    //! The assignment operator
    GlobalIterativeSolverGMRES &operator=( const GlobalIterativeSolverGMRES &other );

    //! The function to set the tolerance
    /*!
     * \param tol Tolerance to be set.
     */
    void setTol( Double tol );

    //! The funciton to set the maximum number of iterations
    /*
     * \para maxIter The maximum number of iterations to be set
     */
    void setMaxIter( UInt maxIter );

    //! The function to set the flag if iteration info should be printed to screen
    /*!
     * \param verbose The flag info to be set
     */
    void setVerbose( bool verbose );

    //! The function to return the residual information
    /*!
     * \param res The vector containing the residual info
     */
    void getLog( vector<Double> &res ) const;

    //! The copy function
    /*!
     * Generates a new copy of this object and returns a pointer to the copy.
     */
    virtual GlobalIterativeSolver *copy() const;

    //! The function to solve the system
    virtual void solve(const vector<ComplexVector> &rhs, vector<ComplexVector> &res ) const;

  private:
    //! An auxilliary function that transforms a vector of vectors over subdomains into a global vector
    /*!
     * \param comm The MPI communicator for the vector storage
     * \param x The input vector to be vectorized
     * \param res The output vector
     *
     * This function merely combines the vectors in x to one big vector res.
     * WARNING: No memory is allocated, memory has to be preallocated!
     */
    void vectorize( MPI_Comm comm, const vector<ComplexVector> &x, ComplexVector &res ) const;
    
    //! An auxilliary function that transforms a global vector to a vector of vectors over subdomains
    /*!
     * \param comm The MPI communicator for the vector storage
     * \param x The input vector to be devectorized
     * \param res The output
     *
     * This function merely splits the vector into vectors defined over subdomains. 
     * WARNING: No memory is allocated, memory has to be preallocated!
     */
    void devectorize( MPI_Comm comm, const ComplexVector &x, vector<ComplexVector> &res ) const;

    //! The auxilliary functiont that realizes an Arnoldi step for the GMRES method
    void arnoldi( MPI_Comm comm, const vector<ComplexMatrix> &Q, UInt k, ComplexVector &h, vector<ComplexVector> &q ) const;
    
    //! The auxilliary functiont that applies a Givens rotation
    void apply_givens_rotation(ComplexVector &h, ComplexVector &cs, ComplexVector &sn, UInt k ) const;

    //! The auxilliary function that computes a Givens rotation
    ComplexVector givens_rotation( Complex v1, Complex v2 ) const;
};

#endif
