#ifndef ASSEMBLERHELMHOLTZ2DFD_H
#define ASSEMBLERHELMHOLTZ2DFD_H

#include<vector>
#include "typedef.h"
#include "AssemblerHelmholtz.h"
#include "Function.h"

//! A class that realizes a finite difference discretization of AssemblerHelmholtz2D
/*!
 * We assume: 
 * - the order of discretization is the same in each direction
 * - the discretization is uniform in each direction
 * - zero Dirichlet boundary conditions are set at the
 *   boundary as dictated by the PMLs
 */
class AssemblerHelmholtz2DFD : public AssemblerHelmholtz
{
  protected:
    //! The interior domain (without the PML region)
    /*!
     * The domain is given in the same way as in Assembler \sa Assembler.
     */
    vector<DoubleVector> mIntDomain;
    
    //! The domain (including the PML region)
    /*!
     * The domain is given in the same way as in Assembler \sa Assembler.
     */
    vector<DoubleVector> mCompDomain;

    //! The number of elements in the x direction 
    UIntVector mN;

    //! The order of the discretization
    UInt mOrder;

  private:
    //! A cache variable containing all DOF points
    mutable vector<DoubleVector> mPts;

  public:
    //! The standard constructor
    /*!
     * Constructs an empty element. The set functions should be used
     * to set the information of the object
     */
    AssemblerHelmholtz2DFD();

    //! The copy constructor
    /*!
     * Properly copies all member variables.
     */
    AssemblerHelmholtz2DFD( const AssemblerHelmholtz2DFD &other );
    
    //! The destructor
    /*!
     * Frees all allocated memory before destruction.
     */
    virtual ~AssemblerHelmholtz2DFD() {}

    //! The copy function
    /*! 
     * Generates a new copy of this object and returns a pointer to the copy.
     */
    virtual Assembler *copy() const;

    //! The assignment operator
    /*!
     * Properly copies \p other into this object
     */
    virtual AssemblerHelmholtz2DFD &operator=( const AssemblerHelmholtz2DFD &other );

    //! The function to set the domain of the assembler
    /*!
     * \sa Assembler
     *
     * The function sets the interior domain. The computational domain is
     * set once the number of elements are set (\sa setN). This defines a 
     * natural order of only setting N after the domain. Therefore, the function
     * clears out any information about N.
     */
    void setDomain( const vector<DoubleVector> &domain );

    //! The function that returns the domain
    /*!
     * \sa Assembler
     */
    vector<DoubleVector> getDomain() const;

    //! The function that sets the number of elements in each direction
    /*!
     * \sa Assembler
     */
    void setN( const UIntVector &N );
    
    //! The function that returns the number of elements in each direction
    /*!
     * \sa Assembler
     */
    UIntVector getN() const;
    
    //! The function that sets the discretization order
    /*!
     * \sa Assembler
     */
    void setOrder( UInt order );

    //! The function that returns the discretization order
    /*!
     * \sa Assembler
     */
    UInt getOrder() const;

    //! The function that returns the mesh-size
    /*!
     * \sa Assembler
     */
    virtual Double getMeshSize() const;
    
    //! The function that returns the dimension 
    /*!
     * \sa Assembler
     */
    virtual UInt getDim() const { return 2; }
    
    //! The function returning the physical points associated with DOFs
    /*!
     * \sa Assembler
     */
    virtual void getDOFPts( vector<DoubleVector> &res ) const;

    //! The function returning the physical points that are associated with DOFs (e.g. boundary points)
    /*!
     * \sa Assembler
     */
    virtual void getNonDOFPhysPts( vector<DoubleVector> &res ) const;

    //! The function returning the indices of the phsical points adjacent to an element
    /*!
     * \sa Assembler
     */
    virtual void getElPhysPtIndices( UInt index, UIntVector &elPts ) const;

    //! The function expanding values associated with DOFs to values defined over all physical points 
    /*!
     * \sa Assembler
     */
    virtual void expandDOFValsToPhysVals( ComplexVector &res ) const;

    //! The function returning the assembled rhs from a \type RHS 
    /*!
     * \sa Assembler
     */
    virtual void assembleRHS( const RHS *rhs, ComplexVector &res ) const;

    //! The function that assembles the system matrix
    /*!
     * \sa Assembler
     */
    virtual void assembleSystemMatrix( ComplexSparseMatrix &res ) const;

    //! The function returning the number of elements in the problem
    /*!
     * \sa Assembler
     */
    virtual UInt getNEls() const;

  private:
    //! A function that returns an approximation of the delta function
    /*!
     * \param x0 The center of the delta function
     * \param x The evaluation point at which the value of the delta
     *          function is returned
     * \return The value of the approximated delta function
     */
    Double deltaFunction( DoubleVector x0, DoubleVector x ) const;
};

#endif
