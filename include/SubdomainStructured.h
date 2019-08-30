#ifndef SUBDOMAINSTRUCTURED_H
#define SUBDOMAINSTRUCTURED_H

#include "Subdomain.h"
#include "Assembler.h"
#include "LocalLinearSolver.h"

class Function;

//! A class that implements a subdomain with an underlying structured grid.
/*!
 * This class assumes that the traces are defined on a structured grid. The traces
 * for higher oder discretizations are defined by taking several rows of DOFs. This is done
 * so that the row of DOFs in trace1 and traceN neighbouring the volume DOFs of this subdomain
 * and the row of DOFs in trace0 and traceNP neighbouring the volume DOFs of the neighbouring
 * subdomain. The DOFs in between trace0 and trace1 (traceN and traceNP) are divided up in
 * two sets of the same size and assigned to trace0 and trace1 (traceN and traceNP) so that
 * each trace has the same number of DOFs.
 *
 * This class also assumes that the domain decomposition is a 2D one.
 */
class SubdomainStructured : public Subdomain
{
  protected:
    bool mHasMatrixCopy;

    //! The index of the subdomain (contains a vector since the indices are given in each direction separately).
    UIntVector mIndex;

    //! The pointer to the assembly routine
    Assembler *mAssembler;

    //! The pointer to the local solver routine
    LocalLinearSolver *mSolver;

    //! The vector that saves the trace type (LeftRight, UpDown,...) for each direction
    vector<string> mTraceTypes;

    //! The vector that saves the index sets divided up in the non-shifted trace and volume DOFs
    vector<vector<UIntVector>> mVolumeIndexSets;

    //! The vector that saves the index sets divided up in the shifted trace and volume DOFs
    vector<vector<UIntVector>> mVolumeIndexSetsShifted;

    //! The vector that saves the element index sets divided up in the non-shifted trace and volume elements
    vector<vector<UIntVector>> mVolumeElIndexSets;
    
    //! The vector that saves the element index sets divided up in the shifted trace and volume elements
    vector<vector<UIntVector>> mVolumeElIndexSetsShifted;

    //! The vector that saves the trace type for the L-shaped traces
    vector<string> mLTraceTypes;
    
    //! The vector that saves the trace indices for the L-shaped traces inherated from both directions non-shifted 
    vector<vector<UIntVector>> mLTraceIndexSets;
    
    //! The vector that saves the trace indices for the L-shaped traces inherated from both directions shifted 
    vector<vector<UIntVector>> mLTraceIndexSetsShifted;
    
    //! The vector that saves the trace indices for the L-shaped traces inherated from the shifted LeftRight direction and the non-shifted UpDown direction 
    vector<vector<UIntVector>> mLTraceIndexSetsShiftedLR;
    
    //! The vector that saves the trace indices for the L-shaped traces inherated from the non-shifted LeftRight direction and the shifted UpDown direction 
    vector<vector<UIntVector>> mLTraceIndexSetsShiftedUD;

  private:
    //! The physical points corresponding to the DOFs
    mutable vector<DoubleVector> mPts;

    //! All physical points in the subdomain
    mutable vector<DoubleVector> mPhysPts;

    //! The cache for the assembled system matrix
    mutable ComplexSparseMatrix *mA;

    //! The cache for the slice of the matrix coming from the bottom trace0 values in the rows and the bottom trace1 values in the columns
    mutable ComplexSparseMatrix mA01_B;

    //! The cache for the slice of the matrix coming from the bottom trace1 values in the rows and the bottom trace0 values in the columns
    mutable ComplexSparseMatrix mA10_B;

    //! The cache for the slice of the matrix coming from the right trace0 values in the rows and the right trace1 values in the columns
    mutable ComplexSparseMatrix mA01_R;

    //! The cache for the slice of the matrix coming from the right trace1 values in the rows and the right trace0 values in the columns
    mutable ComplexSparseMatrix mA10_R;
    
    //! The cache for the slice of the matrix coming from the top trace0 values in the rows and the top trace1 values in the columns
    mutable ComplexSparseMatrix mA01_T;
    
    //! The cache for the slice of the matrix coming from the top trace1 values in the rows and the top trace0 values in the columns
    mutable ComplexSparseMatrix mA10_T;
    
    //! The cache for the slice of the matrix coming from the left trace0 values in the rows and the left trace1 values in the columns
    mutable ComplexSparseMatrix mA01_L;
    
    //! The cache for the slice of the matrix coming from the left trace1 values in the rows and the left trace0 values in the columns
    mutable ComplexSparseMatrix mA10_L;
    
    //! The cache for the slice of the matrix coming from the bottom-left L-shaped trace0 values in the rows and the bottom-left L-shpafted trace1 values in the columns
    mutable ComplexSparseMatrix mA01_BL;
    
    //! The cache for the slice of the matrix coming from the bottom-left L-shaped trace1 values in the rows and the bottom-left L-shpafted trace0 values in the columns
    mutable ComplexSparseMatrix mA10_BL;
    
    //! The cache for the slice of the matrix coming from the bottom-right L-shaped trace0 values in the rows and the bottom-right L-shpafted trace1 values in the columns
    mutable ComplexSparseMatrix mA01_BR;
    
    //! The cache for the slice of the matrix coming from the bottom-right L-shaped trace1 values in the rows and the bottom-right L-shpafted trace0 values in the columns
    mutable ComplexSparseMatrix mA10_BR;
    
    //! The cache for the slice of the matrix coming from the top-right L-shaped trace0 values in the rows and the top-right L-shpafted trace1 values in the columns
    mutable ComplexSparseMatrix mA01_TR;
    
    //! The cache for the slice of the matrix coming from the top-right L-shaped trace1 values in the rows and the top-right L-shpafted trace0 values in the columns
    mutable ComplexSparseMatrix mA10_TR;
    
    //! The cache for the slice of the matrix coming from the top-left L-shaped trace0 values in the rows and the top-left L-shpafted trace1 values in the columns
    mutable ComplexSparseMatrix mA01_TL;
    
    //! The cache for the slice of the matrix coming from the top-left L-shaped trace1 values in the rows and the top-left L-shpafted trace0 values in the columns
    mutable ComplexSparseMatrix mA10_TL;

    //! The cache for the slice of the matrix coming from the shifted bottom trace0 values in the rows and the shifted bottom trace1 values in the columns
    mutable ComplexSparseMatrix mA01Shifted_B;

    //! The cache for the slice of the matrix coming from the shifted bottom trace1 values in the rows and the shifted bottom trace0 values in the columns
    mutable ComplexSparseMatrix mA10Shifted_B;
    
    //! The cache for the slice of the matrix coming from the shifted right trace0 values in the rows and the shifted right trace1 values in the columns
    mutable ComplexSparseMatrix mA01Shifted_R;

    //! The cache for the slice of the matrix coming from the shifted right trace1 values in the rows and the shifted right trace0 values in the columns
    mutable ComplexSparseMatrix mA10Shifted_R;

    //! The cache for the slice of the matrix coming from the shifted right trace1 values in the rows and the shifted right trace0 values in the columns
    mutable ComplexSparseMatrix mA01Shifted_T;

    //! The cache for the slice of the matrix coming from the shifted top trace0 values in the rows and the shifted top trace1 values in the columns
    mutable ComplexSparseMatrix mA10Shifted_T;
    
    //! The cache for the slice of the matrix coming from the shifted top trace1 values in the rows and the shifted top trace0 values in the columns
    mutable ComplexSparseMatrix mA01Shifted_L;

    //! The cache for the slice of the matrix coming from the shifted left trace0 values in the rows and the shifted left trace1 values in the columns
    mutable ComplexSparseMatrix mA10Shifted_L;

    //! The cache for the slice of the matrix coming from the shifted left trace1 values in the rows and the shifted left trace0 values in the columns
    mutable ComplexSparseMatrix mA01Shifted_BL;

    //! The cache for the slice of the matrix coming from the shfited bottom-left L-shaped trace1 values in the rows and the shifted bottom-left L-shpafted trace0 values in the columns
    mutable ComplexSparseMatrix mA10Shifted_BL;
    
    //! The cache for the slice of the matrix coming from the shifted bottom-right L-shaped trace0 values in the rows and the shifted bottom-right L-shpafted trace1 values in the columns
    mutable ComplexSparseMatrix mA01Shifted_BR;
    
    //! The cache for the slice of the matrix coming from the shifted bottom-right L-shaped trace1 values in the rows and the shifted bottom-right L-shpafted trace0 values in the columns
    mutable ComplexSparseMatrix mA10Shifted_BR;
    
    //! The cache for the slice of the matrix coming from the shifted top-right L-shaped trace0 values in the rows and the shifted top-right L-shpafted trace1 values in the columns
    mutable ComplexSparseMatrix mA01Shifted_TR;

    //! The cache for the slice of the matrix coming from the shifted top-right L-shaped trace1 values in the rows and the shifted top-right L-shpafted trace0 values in the columns
    mutable ComplexSparseMatrix mA10Shifted_TR;
    
    //! The cache for the slice of the matrix coming from the shifted top-left L-shaped trace0 values in the rows and the shifted top-left L-shpafted trace1 values in the columns
    mutable ComplexSparseMatrix mA01Shifted_TL;

    //! The cache for the slice of the matrix coming from the shifted top-left L-shaped trace1 values in the rows and the shifted top-left L-shpafted trace0 values in the columns
    mutable ComplexSparseMatrix mA10Shifted_TL;

    //! The cache for the slice of the matrix coming from the bottom-left L-shaped trace0 values shifted in the LeftRight direction in the rows and the bottom-left L-shaped trace1 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA01ShiftedLR_BL;
    
    //! The cache for the slice of the matrix coming from the bottom-left L-shaped trace1 values shifted in the LeftRight direction in the rows and the bottom-left L-shaped trace0 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA10ShiftedLR_BL;

    //! The cache for the slice of the matrix coming from the bottom-right L-shaped trace0 values shifted in the LeftRight direction in the rows and the bottom-right L-shaped trace1 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA01ShiftedLR_BR;
    
    //! The cache for the slice of the matrix coming from the bottom-right L-shaped trace1 values shifted in the LeftRight direction in the rows and the bottom-right L-shaped trace0 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA10ShiftedLR_BR;
    
    //! The cache for the slice of the matrix coming from the top-right L-shaped trace0 values shifted in the LeftRight direction in the rows and the top-right L-shaped trace1 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA01ShiftedLR_TR;
    
    //! The cache for the slice of the matrix coming from the top-right L-shaped trace1 values shifted in the LeftRight direction in the rows and the top-right L-shaped trace0 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA10ShiftedLR_TR;
    
    //! The cache for the slice of the matrix coming from the top-left L-shaped trace0 values shifted in the LeftRight direction in the rows and the top-left L-shaped trace1 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA01ShiftedLR_TL;
    
    //! The cache for the slice of the matrix coming from the top-left L-shaped trace1 values shifted in the LeftRight direction in the rows and the top-left L-shaped trace0 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA10ShiftedLR_TL;
    
    //! The cache for the slice of the matrix coming from the bottom-left L-shaped trace0 values shifted in the LeftRight direction in the rows and the bottom-left L-shaped trace1 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA01ShiftedUD_BL;
    
    //! The cache for the slice of the matrix coming from the bottom-left L-shaped trace1 values shifted in the LeftRight direction in the rows and the bottom-left L-shaped trace0 values shifted in the LeftRight direction in the columns
    mutable ComplexSparseMatrix mA10ShiftedUD_BL;
    
    //! The cache for the slice of the matrix coming from the bottom-right L-shaped trace0 values shifted in the UpDown direction in the rows and the bottom-right L-shaped trace1 values shifted in the UpDown direction in the columns
    mutable ComplexSparseMatrix mA01ShiftedUD_BR;
    
    //! The cache for the slice of the matrix coming from the bottom-right L-shaped trace1 values shifted in the UpDown direction in the rows and the bottom-right L-shaped trace0 values shifted in the UpDown direction in the columns
    mutable ComplexSparseMatrix mA10ShiftedUD_BR;

    //! The cache for the slice of the matrix coming from the top-right L-shaped trace0 values shifted in the UpDown direction in the rows and the top-right L-shaped trace1 values shifted in the UpDown direction in the columns
    mutable ComplexSparseMatrix mA01ShiftedUD_TR;
    
    //! The cache for the slice of the matrix coming from the top-right L-shaped trace1 values shifted in the UpDown direction in the rows and the top-right L-shaped trace0 values shifted in the UpDown direction in the columns
    mutable ComplexSparseMatrix mA10ShiftedUD_TR;
    
    //! The cache for the slice of the matrix coming from the top-left L-shaped trace0 values shifted in the UpDown direction in the rows and the top-left L-shaped trace1 values shifted in the UpDown direction in the columns
    mutable ComplexSparseMatrix mA01ShiftedUD_TL;
    
    //! The cache for the slice of the matrix coming from the top-left L-shaped trace1 values shifted in the UpDown direction in the rows and the top-left L-shaped trace0 values shifted in the UpDown direction in the columns
    mutable ComplexSparseMatrix mA10ShiftedUD_TL;

    //! The cache for the assembled rhs vector using a function 
    mutable ComplexVector mF;

    //! The cache for the mesh-size
    mutable Double mH;

  public:
    //! The standard constructor, constructs an empty object
    SubdomainStructured();

    //! The copy constructor
    SubdomainStructured( const SubdomainStructured &other );

    //! The assignement operator
    SubdomainStructured &operator=( const SubdomainStructured &other );

    //! The destructor
    virtual ~SubdomainStructured();

    //! The copy function
    /*!
     * Generates a new copy of this object and returns a pointer to the copy.
     */
    virtual Subdomain *copy() const;
    
    //! The function to clear the cache
    virtual void clearCache() const;

    //! The function to set the assembly routine, generates a copy of the assembly routine
    void setAssembler( const Assembler *assembler );

    //! The function that sets the interior domain (without the PML)
    void setDomain( const vector<DoubleVector> &domain );
    
    //! The function that returns the interior domain (without the PML)
    vector<DoubleVector> getDomain() const;

    //! The function that sets the number of elements in each direction
    void setN( const UIntVector &n );

    //! The funciotn that returns the number of elements in each direction
    UIntVector getN() const;

    //! The function that sets the order of discretization
    void setOrder( UInt order );

    //! The function that returns the order of discretization
    UInt getOrder() const;

    //! The function to set the local solver routine, generates a copy of the solver routine
    void setLocalLinearSolver( const LocalLinearSolver *solver );
    
    //! The function that sets the locations for the traces
    /*!
     * This function computes the index sets for all trace and volume DOFs. The location of
     * the traces are given in \p lims. Each entry of lims corresponds to one parametric
     * direction in the order LeftRight, UpDown, FrontBack. Each element of \p lims is itself
     * a vector and contains the locations of the traces in the following way:
     * 0: bottom (left,...) boundary of the computational domain
     * 1: trace0 non-shifted
     * 3: trace1 non-shifted
     * 2: trace0 shifted
     * 4: trace1 shifted
     * 5: traceN non-shifted
     * 7: traceNP non-shifted
     * 6: traceN shifted
     * 8: traceNP shifted
     * 9: top (right,...) boundary of the computational domain
     * The locations are given as the location in the parametric direction of the trace. For
     * example, for the LeftRight traces, the corresponding element of \p lims contains the
     * x-coordinate of the trace locations, and for UpDown traces the y-coordinate, etc.
     */
    void setLimits( const vector<DoubleVector> &lims );

    //! The function to set the row index
    virtual void setRowIndex(UInt I);

    //! The function to set the column index
    virtual void setColIndex(UInt J);

    //! The function that returns the row index
    virtual UInt getRowIndex() const;

    //! The function that returns the column index
    virtual UInt getColIndex() const;

    //! The function that returns the number of DOFs in the discretization
    virtual UInt getSize() const;
    
    virtual UInt getDim() const;

    //! The function that returns the physical locations of the DOFs
    virtual void getDOFPts( vector<DoubleVector> &res ) const;

    //! The function that returns all physical points in the subdomain
    /*!
     * The result is ordered to return all DOF points first and then
     * all physical points that are not DOF points.
     *
     * \p res The result
     */
    virtual void getPhysPts( vector<DoubleVector> &res ) const;

    //! The function that returns information of an element
    /*!
     * The function returns a vector of elements. Each element
     * is described by a vector of indices. Each index is a point index
     * pointing to a physical point following the order ing of getPhysPhts.
     *
     * \p elIndex The element index
     * \p res The result
     */
    virtual void getElPhysPtIndices( UInt elIndex, UIntVector &res ) const;

    //! The funciton that expands a vector of values to all physical points
    /*!
     * Following the ordering of getPhysPts, this function just extends
     * the vector by zeros so that it has length equal to the total number
     * of physical points.
     *
     * \p res The result
     */
    virtual void expandDOFValsToPhysVals( ComplexVector &res ) const;

    //! The functiont that assembles the right hand side from an RHS
    /*!
     * \p rhs The RHS
     * \p res The result
     */
    virtual void assembleRHS( const RHS *rhs, ComplexVector &res ) const;

    //! The function that assembles the system matrix and stores it in cache
    virtual void assembleSystemMatrix() const;

    //! The function that assembles the system matrix, stores it in cache, and returns it
    /*!
     * \param res The Matrix where the system matrix is stored, is overwritten
     */
    virtual void assembleSystemMatrix( ComplexSparseMatrix &res ) const;

    //! The function that factorizes the system matrix
    virtual void factorize();

    //! The function that solves the linear system
    /*!
     * For this function, the local linear system first has to be factorized, otherwise this
     * results in an error message.
     * \param x The rhs of the linear system
     * \param res The solution of the system, is overwritten
     */
    virtual void solve( ComplexVector &x, ComplexVector &res ) const;

    //! The function that returns the mesh-size
    virtual Double getMeshSize() const;
    
    virtual UInt getNEls() const;

    //! The function that returns the index sets of the trace and volume DOFs
    /*!
     * This function returns the indices of the DOF sets in the following order:
     * 0: DOFs below (left of,...) trace0
     * 1: trace0 
     * 2: trace1
     * 3: interior volume DOFs
     * 4: traceN
     * 5: traceNP
     * 6: DOFs above (right of,...) traceNP
     *
     * \param traceType The parametric direction of the index sets (LeftRight, UpDown,...)
     * \param shifted The flag if the shifted or non-shifted index sets are asked for
     * \param res The vector that stores the result, is overwritten
     */
    virtual void getIndexSets( string traceType, bool shifted, vector<UIntVector> &res ) const;

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
    virtual void getElIndexSets( string traceType, bool shifted, vector<UIntVector> &res ) const;

    //! The function that applies the local system matrix to a vector
    /*!
     * \param x The vector the system matrix should be applied to
     * \param res The resulting vector (A*x), is overwritten
     */
    virtual void applySystemMatrix( const ComplexVector &x, ComplexVector &res ) const;
 
    //! The function that returns the trace indices of the L traces
    /*!
     * \param traceType The type of the L trace (e.g. BottomLeft, BottomRight, TopRight, TopLeft..)
     * \param shifted The flag if the shifted or non-shfted should be returned
     * \res The resulting index vector (2-element vector, contains the loc0 and loc1 indices for the corresponding trace), is overwritten
     */
    virtual void getLTraceIndexSets( string traceType, bool shiftedLR, bool shiftedUD, vector<UIntVector> &res ) const;

    //! The function that returns the equivalent sources in order to propagate trace info into the subdomain
    /*!
     * \param traceType The type of the L trace (e.g. BottomLeft, BottomRight, TopRight, TopLeft..)
     * \param shifted The flag if the shifted or non-shfted should be used
     * \param vals0 The vector of values for trace0
     * \param vals1 The vector of values for trace1
     * \param res The resulting vector, contains the equivalent sources at the locations of trace0 and trace1, is overwritten
     */
    virtual void getEquivalentSources( string traceType, bool shiftedLR, bool shiftedUD, ComplexVector vals0, ComplexVector vals1, vector<ComplexVector> &res ) const;
    
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
    virtual void extractLTraces( string LTraceType, bool shiftedLR, bool shiftedUD, ComplexVector valsLR0, ComplexVector valsLR1, ComplexVector valsUD0, ComplexVector valsUD1, vector<ComplexVector> &res ) const;

  private:
    //! The function to compute the index sets
    /*!
     * \param traceType The parametric direction of the index sets (LeftRight, UpDown,...)
     * \param lims The locations of the traces, \sa setLimits
     * \return The non-shifted and shifted index sets (in this order).
     */
    vector<vector<UIntVector>> computeVolumeIndexSets(string traceType, const DoubleVector &lims);

    //! The function to compute the element index sets
    /*!
     * \param traceType The parametric direction of the element index sets (LeftRight, UpDown,...)
     * \param volumeInds The volume index sets of the DOFs (obtained from computeVolumeIndexSets)
     */
    vector<UIntVector> computeVolumeElIndexSets( string traceType, const vector<UIntVector> &volumeInds ) const;

    //! An auxilliary function that computes temporary index sets taking the shifted and non-shifted traces into account.
    /*!
     * \param traceType The parametric direction of the index sets (LeftRight, UpDown,...)
     * \param lims The locations of the traces, \sa setLimits
     * \return The auxilliary index sets.
     */
    vector<UIntVector> assignIndexSets(string traceType, const DoubleVector &lims);
   
    //! The function that computes the L trace indices
    /*!
     * \param traceType The type of the L trace (e.g. BottomLeft, BottomRight, TopRight, TopLeft,...)
     * \param shifted The flag if the shifted or non-shifted indices should be computed
     */
    vector<UIntVector> computeLTraceIndexSets(string traceType, bool shiftedLR, bool shiftedUD);
};

#endif
