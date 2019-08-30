#ifndef SUBDOMAIN_H
#define SUBDOMAIN_H

#include "typedef.h"
#include <vector>

using namespace::std;

class LocalLinearSolver;
class RHS;

//! The abstract class for a subdomain
/*!
 * This class manages the discretization within
 * a subdoain and the definition of traces and their
 * extraction. 
 * This class is overloaded by specific subdomain routines. The abstract class is
 * used in other parts of the code (e.g. a domain decomposition).
 */
class Subdomain
{
  public:
    //! The standard constructor, constructs an empty object
    Subdomain() {}
    
    //! The copy constructor
    Subdomain( const Subdomain &other ) {}
    
    //! The destructor
    virtual ~Subdomain() {}

    //! The copy function
    /*!
     * This function copies the (derived) object properly and returns a pointer to
     * the abstract class. This is important to be able to store different types
     * of assembly routines in other parts of the code while not changing the code
     * and consequently makes the code more modular.
     */
    virtual Subdomain *copy() const = 0;
    
    //! The function to clear the cache
    virtual void clearCache() const = 0;
 
    //! The function that sets the limits for the sets of indices determined by the traces
    /*!
     * the limits have to be described in a 10-element vector:
     * Each limit should be in between two rows of degrees of freedom.
     * The limiters are defined to be below (left of) the first row of DOFs
     * corresponding to
     * 0: trace0 non-shifted
     * 1: trace1 non-shifted
     * 2: trace0 shifted
     * 3: trace1 shifted
     * 4: remaining volume DOFs
     * 5: traceN non-shifted
     * 6: traceNP non-shifted
     * 7: traceN shifted
     * 8: traceNP shifted
     * 9: starting row of DOFs above traceNP shifted
     *
     * \p lims The vector of limits
     */
    virtual void setLimits( const vector<DoubleVector> &lims ) = 0;
    
    //! The function that sets the solver used in each cell
    virtual void setLocalLinearSolver( const LocalLinearSolver *solver ) = 0;

    //! The function to set the row index of the subdomain
    virtual void setRowIndex( UInt I ) = 0;
    
    //! The function to set the column index of the subdomain
    virtual void setColIndex( UInt J ) = 0;
    
    //! The function that returns the row index of the subdomain
    virtual UInt getRowIndex() const = 0;
    
    //! The function that returns the column index of the subdomain
    virtual UInt getColIndex() const = 0;
    
    //! The function that returns the size of the subdomain (i.e. total number of DOFs)
    virtual UInt getSize() const = 0;
    
    //! The function that returns the mesh size
    virtual Double getMeshSize() const = 0;

    //! The function that returns the number of elements
    virtual UInt getNEls() const = 0;

    //! The function that returns the dimension
    virtual UInt getDim() const = 0;

    //! The function that returns the physical points associated with the DOFs
    /*!
     * \param res The result vector of physical points
     */
    virtual void getDOFPts( vector<DoubleVector> &res ) const = 0;

    //! The function that returns all physical points in the subdomain
    /*!
     * The result is ordered to return all DOF points first and then
     * all physical points that are not DOF points.
     *
     * \p res The result
     */
    virtual void getPhysPts( vector<DoubleVector> &res ) const = 0;

    //! The function that returns information of an element
    /*!
     * The function returns a vector of elements. Each element
     * is described by a vector of indices. Each index is a point index
     * pointing to a physical point following the order ing of getPhysPhts.
     *
     * \p elIndex The element index
     * \p res The result
     */
    virtual void getElPhysPtIndices( UInt elIndex, UIntVector &res ) const = 0;

    //! The funciton that expands a vector of values to all physical points
    /*!
     * Following the ordering of getPhysPts, this function just extends
     * the vector by zeros so that it has length equal to the total number
     * of physical points.
     *
     * \p res The result
     */
    virtual void expandDOFValsToPhysVals( ComplexVector &res ) const = 0;

    //! The functiont that assembles the right hand side from an RHS
    /*!
     * \p rhs The RHS
     * \p res The result
     */
    virtual void assembleRHS( const RHS *rhs, ComplexVector &res ) const = 0;
    
    //! The Function that assembles the system matrix,the system matrix is
    //  stored in cache.
    /*!
     * \param res The resulting system matrix, is overwritten
     */
    virtual void assembleSystemMatrix( ComplexSparseMatrix &res ) const = 0;
    
    //! The Function that assembles the system matrix and stores it in cache.
    virtual void assembleSystemMatrix() const = 0;

    //! The Function that factorizes the system matrix in the subdomain
    virtual void factorize() = 0;
    
    //! The function that solves the linear system in the subdomain
    /*!
     * \param x The right hand side for the linear system
     * \param res The result of the solved system, is overwritten
     */
    virtual void solve( ComplexVector &x, ComplexVector &res ) const = 0;
    
    //! The function that returns the index sets
    /*!
     * Returns the index sets after ordering them corresponding to volume indices
     * and traces. The order is:
     * 0: bottom DOFs (below trace0)
     * 1: trace0
     * 2: trace1
     * 3: interior DOFS (between trace1 and traceN)
     * 4: traceN
     * 5: traceNP
     * 6: top DOFs (above traceNP)
     * \param traceType The type of the trace (e.g. UpDown, LeftRight...)
     * \param shifted Flag if the non-shifted or shifted indices should be used
     * \param res The resulting sets of indices
     */
    virtual void getIndexSets( string traceType, bool shifted, vector<UIntVector> &res ) const = 0;

    //! The function that returns the element index sets
    /*!
     * Returns the index sets of elements similar to the DOF points.
     * The order is:
     * 0: All elements below trace0
     * 1: All elements in trace0 (might be empty, depending on the trace thickness)
     * 2: All elements between trace0 and trace1
     * 3: All elements in trace1 (might be empty, depending on the trace thickness)
     * 4: All elements between trace1 and traceN
     * 5: All elements in traceN (might be empty, depending on the trace thickness)
     * 6: All elements between traceN and traceNP
     * 7: All elements in traceNP (might be empty, depending on the trace thickness)
     * 8: All elements above traceNP
     * \p traceType The type of the trace (e.g. UpDown, LeftRight,...)
     * \p shifted The flag if the non-shifted or shifted traces should be used
     * \p res The result
     */
    virtual void getElIndexSets( string traceType, bool shifted, vector<UIntVector> &res ) const = 0;

    //! The function that applies the system matrix
    /*!
     * \param x The vector the system matrix should be applied to
     * \param res The resulting vector (res=A*x), is overwritten
     */
    virtual void applySystemMatrix( const ComplexVector &x, ComplexVector &res ) const = 0;
   
    //! The function that returns the trace indices of the L traces
    /*!
     * \param traceType The type of the L trace (e.g. BottomLeft, BottomRight, TopRight, TopLeft..)
     * \param shifted The flag if the shifted or non-shfted should be returned
     * \res The resulting index vector (2-element vector, contains the loc0 and loc1 indices for the corresponding trace), is overwritten
     */
    virtual void getLTraceIndexSets( string traceType, bool shiftedLR, bool shiftedUD, vector<UIntVector> &res ) const = 0;

    //! The function that returns the equivalent sources in order to propagate trace info into the subdomain
    /*!
     * \param traceType The type of the L trace (e.g. BottomLeft, BottomRight, TopRight, TopLeft..)
     * \param shifted The flag if the shifted or non-shfted should be used
     * \param vals0 The vector of values for trace0
     * \param vals1 The vector of values for trace1
     * \param res The resulting vector, contains the equivalent sources at the locations of trace0 and trace1, is overwritten
     */
    virtual void getEquivalentSources( string traceType, bool shiftedLR, bool shiftedUD, ComplexVector vals0, ComplexVector vals1, vector<ComplexVector> &res ) const = 0;

    //! The function that extracts L trace values from standard (horizontal or vertical) traces
    /*!
     * \param traceType The type of the L trace (e.g. BottomLeft, BottomRight, TopRight, TopLeft..)
     * \param shifted The flag if the shifted or non-shfted should be used
     * \param valsLR0 The trace0 values of the LeftRight trace
     * \param valsLR1 The trace1 values of the LeftRight trace
     * \param valsUD0 The trace0 values of the UpDown trace
     * \param valsUD1 The trace1 values of the UpDown trace
     * \param res The resulting vector, contains the L trace values of trace0 and trace1, is overwritten
     */
    virtual void extractLTraces( string LTraceType, bool shiftedLR, bool shiftedUD, ComplexVector valsLR0, ComplexVector valsLR1, ComplexVector valsUD0, ComplexVector valsUD1, vector<ComplexVector> &res ) const = 0;
};

#endif
