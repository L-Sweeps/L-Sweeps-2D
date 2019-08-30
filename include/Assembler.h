#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "typedef.h"

class RHS;

//! The abstract class for an assembly routine
/*!
 * This class is overloaded by specific assembly routines.
 */
class Assembler
{
  public:
    //! The standard constructor
    Assembler() {}

    //! The copy constructor
    Assembler( const Assembler &other ) {}

    //! The destructor
    virtual ~Assembler() {}


    //! The copy function
    /*!
     * This function copies the (derived) object properly and returns a pointer to
     * the abstract class. This is important to be able to store different types
     * of assembly routines in other parts of the code while not changing the code
     * and consequently makes the code more modular.
     */
    virtual Assembler *copy() const = 0;

    //! A function to set the domain of the assmebler
    /*!
     * \p dom The domain info that should be used, the domain is assumed to be square
     *
     * The domain is given as a vector of vector of doubles. The vector has length equal to the dimension.
     * Each element (vector of doubles) contains two values determining the interval of the domain in this direction.
     */
    virtual void setDomain(const vector<DoubleVector> &dom) = 0;

    //! A function returning the domain
    virtual vector<DoubleVector> getDomain() const = 0;

    //! A function that sets the number of elements used for discretization
    /*!
     * \p n The vector holding the number of elements
     *
     * The elements are given by a vector of length equal to the dimension. Each element of the vector
     * determines the number of elements used in this direction.
     */
    virtual void setN( const UIntVector &n ) = 0;
    
    //! A function returning the number of elements in each direction as a vector
    virtual UIntVector getN() const = 0;

    //! A function to set the order of the discretization
    /*!
     * \p order The order to be set given as an integer
     */
    virtual void setOrder( UInt order ) = 0;
    
    //! A function returning the order
    virtual UInt getOrder() const = 0;

    //! A function returning the mesh size of the underlying discretization 
    virtual Double getMeshSize() const = 0;

    //! A function returning the dimension of the underlying discretized problem
    virtual UInt getDim() const = 0;

    //! A function returning the physical points associated with DOFs
    /*!
     * \p res is overwritten with the returned physical ponts.
     */
    virtual void getDOFPts( vector<DoubleVector> &res ) const = 0;

    //! A function returning the physical points that are associated with DOFs (e.g. boundary points)
    /*!
     * \p res is overwritten with the returned physical ponts.
     */
    virtual void getNonDOFPhysPts( vector<DoubleVector> &res ) const = 0;

    //! A function returning the indices of the phsical points adjacent to an element
    /*!
     * The physical points pointed to include DOF and non-DOF points.
     *
     * \p index The index of the element
     * \p elPts is overwritten with the returned indices to physical points
     */
    virtual void getElPhysPtIndices( UInt index, UIntVector &elPts ) const = 0;

    //! A function expanding values associated with DOFs to values defined over all physical points 
    /*!
     * This function takes a vector of values associated with the DOFs of the problem.
     * The input vector is ordered so that it is consistent with getDOFPts().
     * The values are expanded by zeros to include all physical points that are not DOFs.
     * All added values are appended at the end of the input vector.
     *
     * \p res The input and output vector.
     */
    virtual void expandDOFValsToPhysVals( ComplexVector &res ) const = 0;


    //! A function returning the assembled rhs from a \type RHS 
    /*!
     * \param rhs The rhight hand side used to assemble the rhs.
     * \param res The vector overwritten with the assembled rhs.
     */

    virtual void assembleRHS( const RHS *rhs, ComplexVector &res ) const = 0;

    //! A function returning the assembled system matrix
    /*!
     * \param res The matrix overwritten with the assembled system matrix
     */
    virtual void assembleSystemMatrix( ComplexSparseMatrix &res ) const = 0;

    //! A function returning the number of elements in the problem
    virtual UInt getNEls() const = 0;
};


#endif

