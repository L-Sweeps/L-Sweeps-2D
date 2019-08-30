#include "SubdomainStructured.h"
#include "PMLCubic2D.h"
#include "AssemblerHelmholtz2DFD.h"
#include <iostream>
#include "utils.h"

using namespace::std;
using namespace::Eigen;

// The standard constructor generating an empty object
// We make sure that NULL vectors are properly generated and the index vector has the 
// correct size
SubdomainStructured::SubdomainStructured() : mAssembler(NULL), mSolver(NULL), mIndex(2), mPhysPts(0), mPts(0), mA(NULL), mHasMatrixCopy(true)
{

}

// The destructor
// The assembler and local solver routines are freed because they are copied in the process.
SubdomainStructured::~SubdomainStructured()
{
  if( mAssembler!=NULL )
    delete mAssembler;

  if( mSolver!=NULL )
    delete mSolver;

  if( mA!=NULL && mHasMatrixCopy )
    delete mA;
}

// The copy constructor
SubdomainStructured::SubdomainStructured( const SubdomainStructured &other ): mHasMatrixCopy(other.mHasMatrixCopy)
{
  // copy the member variables
  mIndex = other.mIndex;
  if( other.mAssembler != NULL )
    mAssembler = other.mAssembler->copy();
  else
    mAssembler = NULL;
  if( other.mSolver != NULL )
    mSolver = other.mSolver->copy();
  else
    mSolver = NULL;

  mA = other.mA;
  if( other.mA != NULL && mHasMatrixCopy )
    mA = new ComplexSparseMatrix(*other.mA);

  mTraceTypes = other.mTraceTypes;
  mLTraceTypes = other.mLTraceTypes;
  mVolumeIndexSets = other.mVolumeIndexSets;
  mVolumeElIndexSets = other.mVolumeElIndexSets;
  mVolumeIndexSetsShifted = other.mVolumeIndexSetsShifted;
  mVolumeElIndexSetsShifted = other.mVolumeElIndexSetsShifted;
  mLTraceIndexSets = other.mLTraceIndexSets;
  mLTraceIndexSetsShifted = other.mLTraceIndexSetsShifted;
  mLTraceIndexSetsShiftedLR = other.mLTraceIndexSetsShiftedLR;
  mLTraceIndexSetsShiftedUD = other.mLTraceIndexSetsShiftedUD;

  // copy the cache
  mPts = other.mPts;
  mPhysPts = other.mPhysPts;
  mA01_B = other.mA01_B;
  mA10_B = other.mA10_B;
  mA01_R = other.mA01_R;
  mA10_R = other.mA10_R;
  mA01_T = other.mA01_T;
  mA10_T = other.mA10_T;
  mA01_L = other.mA01_L;
  mA10_L = other.mA10_L;
  mA01_BL = other.mA01_BL;
  mA10_BL = other.mA10_BL;
  mA01_BR = other.mA01_BR;
  mA10_BR = other.mA10_BR;
  mA01_TR = other.mA01_TR;
  mA10_TR = other.mA10_TR;
  mA01_TL = other.mA01_TL;
  mA10_TL = other.mA10_TL;
  mA01Shifted_B = other.mA01Shifted_B;
  mA10Shifted_B = other.mA10Shifted_B;
  mA01Shifted_R = other.mA01Shifted_R;
  mA10Shifted_R = other.mA10Shifted_R;
  mA01Shifted_T = other.mA01Shifted_T;
  mA10Shifted_T = other.mA10Shifted_T;
  mA01Shifted_L = other.mA01Shifted_L;
  mA10Shifted_L = other.mA10Shifted_L;
  mA01Shifted_BL = other.mA01Shifted_BL;
  mA10Shifted_BL = other.mA10Shifted_BL;
  mA01Shifted_BR = other.mA01Shifted_BR;
  mA10Shifted_BR = other.mA10Shifted_BR;
  mA01Shifted_TR = other.mA01Shifted_TR;
  mA10Shifted_TR = other.mA10Shifted_TR;
  mA01Shifted_TL = other.mA01Shifted_TL;
  mA10Shifted_TL = other.mA10Shifted_TL;
  mA01ShiftedLR_BL = other.mA01ShiftedLR_BL;
  mA10ShiftedLR_BL = other.mA10ShiftedLR_BL;
  mA01ShiftedLR_BR = other.mA01ShiftedLR_BR;
  mA10ShiftedLR_BR = other.mA10ShiftedLR_BR;
  mA01ShiftedLR_TR = other.mA01ShiftedLR_TR;
  mA10ShiftedLR_TR = other.mA10ShiftedLR_TR;
  mA01ShiftedLR_TL = other.mA01ShiftedLR_TL;
  mA10ShiftedLR_TL = other.mA10ShiftedLR_TL;
  mA01ShiftedUD_BL = other.mA01ShiftedUD_BL;
  mA10ShiftedUD_BL = other.mA10ShiftedUD_BL;
  mA01ShiftedUD_BR = other.mA01ShiftedUD_BR;
  mA10ShiftedUD_BR = other.mA10ShiftedUD_BR;
  mA01ShiftedUD_TR = other.mA01ShiftedUD_TR;
  mA10ShiftedUD_TR = other.mA10ShiftedUD_TR;
  mA01ShiftedUD_TL = other.mA01ShiftedUD_TL;
  mA10ShiftedUD_TL = other.mA10ShiftedUD_TL;
  mF = other.mF;
  mH = other.mH;

  // if there is a solver set in other, we have to assign the matrix
  // in order for it to work (Eigen does not allow for copy constructors,
  // that's why we have to do this for general solvers too in order to keep
  // the code consistent).
  if( mSolver != NULL && mA!=NULL )
  {
    mSolver->setMatrix(mA);
  }
}

// The assignment operator
SubdomainStructured &SubdomainStructured::operator=( const SubdomainStructured &other )
{
  // only do soething if other is not the same as this object
  if( this != &other )
  {
    mHasMatrixCopy = other.mHasMatrixCopy;

    // if assembler and local solver routines have been previously set, free the memory
    if( mAssembler != NULL )
      delete mAssembler;
    if( mSolver != NULL )
      delete mSolver;
    if( mA != NULL && mHasMatrixCopy )
      delete mA;

    // assign the member variables (copy assembler and solver)
    mIndex = other.mIndex;
    if( other.mAssembler != NULL )
      mAssembler = other.mAssembler->copy();
    else
      mAssembler = NULL;
    if( other.mSolver != NULL )
      mSolver = other.mSolver->copy();
    else
      mSolver = NULL;
    mA = other.mA;
    if( other.mA != NULL )
      mA = new ComplexSparseMatrix(*other.mA);

    mTraceTypes = other.mTraceTypes;
    mLTraceTypes = other.mLTraceTypes;
    mVolumeIndexSets = other.mVolumeIndexSets;
    mVolumeElIndexSets = other.mVolumeElIndexSets;
    mVolumeIndexSetsShifted = other.mVolumeIndexSetsShifted;
    mVolumeElIndexSetsShifted = other.mVolumeElIndexSetsShifted;
    mLTraceIndexSets = other.mLTraceIndexSets;
    mLTraceIndexSetsShifted = other.mLTraceIndexSetsShifted;
    mLTraceIndexSetsShiftedLR = other.mLTraceIndexSetsShiftedLR;
    mLTraceIndexSetsShiftedUD = other.mLTraceIndexSetsShiftedUD;

    // assign the cache
    mPts = other.mPts;
    mPhysPts = other.mPhysPts;
    mA01_B = other.mA01_B;
    mA10_B = other.mA10_B;
    mA01_R = other.mA01_R;
    mA10_R = other.mA10_R;
    mA01_T = other.mA01_T;
    mA10_T = other.mA10_T;
    mA01_L = other.mA01_L;
    mA10_L = other.mA10_L;
    mA01_BL = other.mA01_BL;
    mA10_BL = other.mA10_BL;
    mA01_BR = other.mA01_BR;
    mA10_BR = other.mA10_BR;
    mA01_TR = other.mA01_TR;
    mA10_TR = other.mA10_TR;
    mA01_TL = other.mA01_TL;
    mA10_TL = other.mA10_TL;
    mA01Shifted_B = other.mA01Shifted_B;
    mA10Shifted_B = other.mA10Shifted_B;
    mA01Shifted_R = other.mA01Shifted_R;
    mA10Shifted_R = other.mA10Shifted_R;
    mA01Shifted_T = other.mA01Shifted_T;
    mA10Shifted_T = other.mA10Shifted_T;
    mA01Shifted_L = other.mA01Shifted_L;
    mA10Shifted_L = other.mA10Shifted_L;
    mA01Shifted_BL = other.mA01Shifted_BL;
    mA10Shifted_BL = other.mA10Shifted_BL;
    mA01Shifted_BR = other.mA01Shifted_BR;
    mA10Shifted_BR = other.mA10Shifted_BR;
    mA01Shifted_TR = other.mA01Shifted_TR;
    mA10Shifted_TR = other.mA10Shifted_TR;
    mA01Shifted_TL = other.mA01Shifted_TL;
    mA10Shifted_TL = other.mA10Shifted_TL;
    mA01ShiftedLR_BL = other.mA01ShiftedLR_BL;
    mA10ShiftedLR_BL = other.mA10ShiftedLR_BL;
    mA01ShiftedLR_BR = other.mA01ShiftedLR_BR;
    mA10ShiftedLR_BR = other.mA10ShiftedLR_BR;
    mA01ShiftedLR_TR = other.mA01ShiftedLR_TR;
    mA10ShiftedLR_TR = other.mA10ShiftedLR_TR;
    mA01ShiftedLR_TL = other.mA01ShiftedLR_TL;
    mA10ShiftedLR_TL = other.mA10ShiftedLR_TL;
    mA01ShiftedUD_BL = other.mA01ShiftedUD_BL;
    mA10ShiftedUD_BL = other.mA10ShiftedUD_BL;
    mA01ShiftedUD_BR = other.mA01ShiftedUD_BR;
    mA10ShiftedUD_BR = other.mA10ShiftedUD_BR;
    mA01ShiftedUD_TR = other.mA01ShiftedUD_TR;
    mA10ShiftedUD_TR = other.mA10ShiftedUD_TR;
    mA01ShiftedUD_TL = other.mA01ShiftedUD_TL;
    mA10ShiftedUD_TL = other.mA10ShiftedUD_TL;
    mF = other.mF;
    mH = other.mH;

    // if there is a solver set in other, we have to assign the matrix
    // in order for it to work (Eigen does not allow for copy constructors,
    // that's why we have to do this for general solvers too in order to keep
    // the code consistent).
    if( mSolver != NULL && mA != NULL )
    {
      mSolver->setMatrix(mA);
    }
  }
  return *this;
}

// the copy function, generates a copy and returns a pointer to it
Subdomain *SubdomainStructured::copy() const
{
  return new SubdomainStructured(*this);
}

// clears the cache
void SubdomainStructured::clearCache() const
{
  if( mA!=NULL && mHasMatrixCopy )
  {
    delete mA;
  }
  mA=NULL;
  mPts.clear();
  mPhysPts.clear();
  mA01_B.resize(0,0);
  mA10_B.resize(0,0);
  mA01_R.resize(0,0);
  mA10_R.resize(0,0);
  mA01_T.resize(0,0);
  mA10_T.resize(0,0);
  mA01_L.resize(0,0);
  mA10_L.resize(0,0);
  mA01_BL.resize(0,0);
  mA10_BL.resize(0,0);
  mA01_BR.resize(0,0);
  mA10_BR.resize(0,0);
  mA01_TR.resize(0,0);
  mA10_TR.resize(0,0);
  mA01_TL.resize(0,0);
  mA10_TL.resize(0,0);
  mA01Shifted_B.resize(0,0);
  mA10Shifted_B.resize(0,0);
  mA01Shifted_R.resize(0,0);
  mA10Shifted_R.resize(0,0);
  mA01Shifted_T.resize(0,0);
  mA10Shifted_T.resize(0,0);
  mA01Shifted_L.resize(0,0);
  mA10Shifted_L.resize(0,0);
  mA01Shifted_BL.resize(0,0);
  mA10Shifted_BL.resize(0,0);
  mA01Shifted_BR.resize(0,0);
  mA10Shifted_BR.resize(0,0);
  mA01Shifted_TR.resize(0,0);
  mA10Shifted_TR.resize(0,0);
  mA01Shifted_TL.resize(0,0);
  mA10Shifted_TL.resize(0,0);
  mA01ShiftedLR_BL.resize(0,0);
  mA10ShiftedLR_BL.resize(0,0);
  mA01ShiftedLR_BR.resize(0,0);
  mA10ShiftedLR_BR.resize(0,0);
  mA01ShiftedLR_TR.resize(0,0);
  mA10ShiftedLR_TR.resize(0,0);
  mA01ShiftedLR_TL.resize(0,0);
  mA10ShiftedLR_TL.resize(0,0);
  mA01ShiftedUD_BL.resize(0,0);
  mA10ShiftedUD_BL.resize(0,0);
  mA01ShiftedUD_BR.resize(0,0);
  mA10ShiftedUD_BR.resize(0,0);
  mA01ShiftedUD_TR.resize(0,0);
  mA10ShiftedUD_TR.resize(0,0);
  mA01ShiftedUD_TL.resize(0,0);
  mA10ShiftedUD_TL.resize(0,0);
  mF.resize(0);
  mH = 0.0;
}

// sets the assembler routine, the cache is cleared in the process
void SubdomainStructured::setAssembler( const Assembler *assembler )
{
  // clear cache
  clearCache();

  // if there was an assembler routine previously set, free the memory
  if( mAssembler != NULL )
    delete mAssembler;

  // copy the assembler routine
  mAssembler = assembler->copy();
}

// sets the limits of the traces and computes the trace and volume index sets
void SubdomainStructured::setLimits( const vector<DoubleVector> &lims )
{
  // get the dimension of the discretization
  UInt dim = mAssembler->getDim();
  // traces have to be defined in each parametric direction 
  assert(lims.size() == dim);

  // the limits have to be described in a 10-element vector:
  // Each limit should be in between two rows of degrees of freedom.
  // The limiters are defined to be below (left of) the first row of DOFs
  // corresponding to
  // 0: trace0 non-shifted
  // 1: trace1 non-shifted
  // 2: trace0 shifted
  // 3: trace1 shifted
  // 4: remaining volume DOFs
  // 5: traceN non-shifted
  // 6: traceNP non-shifted
  // 7: traceN shifted
  // 8: traceNP shifted
  // 9: starting row of DOFs above traceNP shifted
  for( UInt k=0; k<dim; k++ )
    assert(lims[k].rows() == 10 );

  // clear cache
  clearCache();

  // so far we only support 2D domain decompositions
  // set the trace types, compute the index sets and assign them
  mTraceTypes.push_back("LeftRight");
  vector<vector<UIntVector>> temp = computeVolumeIndexSets("LeftRight",lims[0]);
  mVolumeIndexSets.push_back( temp[0] );
  mVolumeIndexSetsShifted.push_back( temp[1] );
  mVolumeElIndexSets.push_back( computeVolumeElIndexSets("LeftRight", temp[0] ) ); 
  mVolumeElIndexSetsShifted.push_back( computeVolumeElIndexSets("LeftRight", temp[1] ) ); 

  mTraceTypes.push_back("UpDown");
  temp = computeVolumeIndexSets("UpDown",lims[1]);
  mVolumeIndexSets.push_back( temp[0] );
  mVolumeIndexSetsShifted.push_back( temp[1] );
  mVolumeElIndexSets.push_back( computeVolumeElIndexSets("UpDown", temp[0] ) ); 
  mVolumeElIndexSetsShifted.push_back( computeVolumeElIndexSets("UpDown", temp[1] ) ); 

  mLTraceTypes.push_back("BottomLeft");
  temp[0] = computeLTraceIndexSets("BottomLeft", false, false);
  mLTraceIndexSets.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("BottomLeft", true, true);
  mLTraceIndexSetsShifted.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("BottomLeft", true, false);
  mLTraceIndexSetsShiftedLR.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("BottomLeft", false, true);
  mLTraceIndexSetsShiftedUD.push_back( temp[0] );

  mLTraceTypes.push_back("BottomRight");
  temp[0] = computeLTraceIndexSets("BottomRight", false, false);
  mLTraceIndexSets.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("BottomRight", true, true);
  mLTraceIndexSetsShifted.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("BottomRight", true, false);
  mLTraceIndexSetsShiftedLR.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("BottomRight", false, true);
  mLTraceIndexSetsShiftedUD.push_back( temp[0] );

  mLTraceTypes.push_back("TopRight");
  temp[0] = computeLTraceIndexSets("TopRight", false, false);
  mLTraceIndexSets.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("TopRight", true, true);
  mLTraceIndexSetsShifted.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("TopRight", true, false);
  mLTraceIndexSetsShiftedLR.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("TopRight", false, true);
  mLTraceIndexSetsShiftedUD.push_back( temp[0] );

  mLTraceTypes.push_back("TopLeft");
  temp[0] = computeLTraceIndexSets("TopLeft", false, false);
  mLTraceIndexSets.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("TopLeft", true, true);
  mLTraceIndexSetsShifted.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("TopLeft", true, false);
  mLTraceIndexSetsShiftedLR.push_back( temp[0] );
  temp[0] = computeLTraceIndexSets("TopLeft", false, true);
  mLTraceIndexSetsShiftedUD.push_back( temp[0] );
}

// sets the local solver, clears cache
void SubdomainStructured::setLocalLinearSolver( const LocalLinearSolver *solver )
{
  // clear cache
  clearCache();

  // if a solver had been previously set, free the memory
  if( mSolver != NULL )
    delete mSolver;

  // copy the solver
  mSolver = solver->copy();

  // if a system matrix had been previously set, set the matrix for the solver
  if( mA!=NULL )
  {
    mSolver->setMatrix(mA);
  }
}

// sets the domain
void SubdomainStructured::setDomain( const vector<DoubleVector> &domain )
{
  mAssembler->setDomain(domain);
}

// returns the domain
vector<DoubleVector> SubdomainStructured::getDomain() const
{
  return mAssembler->getDomain();
}

// sets the number of elements used for discretization in each direction
void SubdomainStructured::setN( const UIntVector &n )
{
  mAssembler->setN(n);
}

// returns the number of elements used for discretization in each direction
UIntVector SubdomainStructured::getN() const
{
  return mAssembler->getN();
}

// sets the order of discretization
void SubdomainStructured::setOrder( UInt order )
{
  mAssembler->setOrder(order);
}

// returns the order of discretization
UInt SubdomainStructured::getOrder() const
{
  return mAssembler->getOrder();
}

// sets the row index
void SubdomainStructured::setRowIndex( UInt I )
{
  mIndex[0] = I;
}

// sets the column index
void SubdomainStructured::setColIndex( UInt J )
{
  mIndex[1] = J;
}

// returns the row index
UInt SubdomainStructured::getRowIndex() const
{
  return mIndex[0];
}

// returns the column index
UInt SubdomainStructured::getColIndex() const
{
  return mIndex[1];
}

// returns the number of DOFs
UInt SubdomainStructured::getSize() const
{
  // if the physical DOF points have not been computed yet, compute them
  if( mPts.size()==0 )
    mAssembler->getDOFPts( mPts );

  // return the number of DOFs
  return mPts.size();
}

// returns the dimension of the problem
UInt SubdomainStructured::getDim() const
{
  return mAssembler->getDim();
}

// returns the physical points corresponding to the DOFs
void SubdomainStructured::getDOFPts( vector<DoubleVector> &res ) const
{
  // if the physical DOF points have not been computed yet, compute them
  if( mPts.size()==0 )
    mAssembler->getDOFPts( mPts );

  // copy the points into the result
  res = mPts;
}

// returns all physical points in the problem
void SubdomainStructured::getPhysPts( vector<DoubleVector> &res ) const
{
  // if the physical points are zero, compute them
  if( mPhysPts.size()==0 )
  {
    // get the DOF points
    getDOFPts( mPhysPts );
    // get the non-DOF points
    vector<DoubleVector> temp;
    mAssembler->getNonDOFPhysPts( temp );
    // combine them
    mPhysPts.insert(mPhysPts.end(), temp.begin(), temp.end());
  }
  // copy the physical points into the result
  res = mPhysPts;
}

// get the element to point information, uses the same routine from the assembly routine
void SubdomainStructured::getElPhysPtIndices( UInt index, UIntVector &elPts ) const
{
  return mAssembler->getElPhysPtIndices( index, elPts );
}

// expands values defined on the DOF points to all physical points, simply uses the equivalent routine from the assembler
void SubdomainStructured::expandDOFValsToPhysVals( ComplexVector &res ) const
{
  mAssembler->expandDOFValsToPhysVals( res ); 
}

// assembles the right hand side from a given RHS
void SubdomainStructured::assembleRHS( const RHS *rhs, ComplexVector &res ) const
{
  // if no rhs is in the cached yet, compute it
  if( mF.rows()==0 )
    mAssembler->assembleRHS(rhs,mF);

  // return the rhs
  res = mF;
}

// assembles the system atrix and stores it in cahce
void SubdomainStructured::assembleSystemMatrix() const
{
  if( !mHasMatrixCopy )
  {
    cerr << "ERROR: Matrix can only be assembled for subdomains that store system matrices!" << endl;
    exit(1);
  }
  // if matrix is not assembled yet, assemble it and cache it
  if( mA==NULL )
  {
    ComplexSparseMatrix A;
    mAssembler->assembleSystemMatrix(A);
    mA = new ComplexSparseMatrix(A);
  }

  // compress the matrix
  mA->makeCompressed();

  // get the trace indices and set up the matrices needed for the equivalent sources
  vector<UIntVector> indexSets1, indexSets2;
  vector<UIntVector> indexSets1Shifted, indexSets2Shifted;
  getIndexSets("LeftRight",false,indexSets1);
  getIndexSets("UpDown",false,indexSets2);
  getIndexSets("LeftRight",true,indexSets1Shifted);
  getIndexSets("UpDown",true,indexSets2Shifted);

  vector<UIntVector> LTraceIndexSetsBL, LTraceIndexSetsBR, LTraceIndexSetsTR, LTraceIndexSetsTL;
  vector<UIntVector> LTraceIndexSetsBLShifted, LTraceIndexSetsBRShifted, LTraceIndexSetsTRShifted, LTraceIndexSetsTLShifted;
  vector<UIntVector> LTraceIndexSetsBLShiftedLR, LTraceIndexSetsBRShiftedLR, LTraceIndexSetsTRShiftedLR, LTraceIndexSetsTLShiftedLR;
  vector<UIntVector> LTraceIndexSetsBLShiftedUD, LTraceIndexSetsBRShiftedUD, LTraceIndexSetsTRShiftedUD, LTraceIndexSetsTLShiftedUD;
  getLTraceIndexSets("BottomLeft",false,false,LTraceIndexSetsBL);
  getLTraceIndexSets("BottomRight",false,false,LTraceIndexSetsBR);
  getLTraceIndexSets("TopRight",false,false,LTraceIndexSetsTR);
  getLTraceIndexSets("TopLeft",false,false,LTraceIndexSetsTL);
  getLTraceIndexSets("BottomLeft",true,true,LTraceIndexSetsBLShifted);
  getLTraceIndexSets("BottomRight",true,true,LTraceIndexSetsBRShifted);
  getLTraceIndexSets("TopRight",true,true,LTraceIndexSetsTRShifted);
  getLTraceIndexSets("TopLeft",true,true,LTraceIndexSetsTLShifted);
  getLTraceIndexSets("BottomLeft",true,false,LTraceIndexSetsBLShiftedLR);
  getLTraceIndexSets("BottomRight",true,false,LTraceIndexSetsBRShiftedLR);
  getLTraceIndexSets("TopRight",true,false,LTraceIndexSetsTRShiftedLR);
  getLTraceIndexSets("TopLeft",true,false,LTraceIndexSetsTLShiftedLR);
  getLTraceIndexSets("BottomLeft",false,true,LTraceIndexSetsBLShiftedUD);
  getLTraceIndexSets("BottomRight",false,true,LTraceIndexSetsBRShiftedUD);
  getLTraceIndexSets("TopRight",false,true,LTraceIndexSetsTRShiftedUD);
  getLTraceIndexSets("TopLeft",false,true,LTraceIndexSetsTLShiftedUD);

  mA01_B = slice(*mA,indexSets2[1],indexSets2[2]); 
  mA10_B = slice(*mA,indexSets2[2],indexSets2[1]); 
  mA01_R = slice(*mA,indexSets1[5],indexSets1[4]); 
  mA10_R = slice(*mA,indexSets1[4],indexSets1[5]); 
  mA01_T = slice(*mA,indexSets2[5],indexSets2[4]); 
  mA10_T = slice(*mA,indexSets2[4],indexSets2[5]); 
  mA01_L = slice(*mA,indexSets1[1],indexSets1[2]); 
  mA10_L = slice(*mA,indexSets1[2],indexSets1[1]); 
  mA01_BL = slice(*mA,LTraceIndexSetsBL[0],LTraceIndexSetsBL[1]);
  mA10_BL = slice(*mA,LTraceIndexSetsBL[1],LTraceIndexSetsBL[0]);
  mA01_BR = slice(*mA,LTraceIndexSetsBR[0],LTraceIndexSetsBR[1]);
  mA10_BR = slice(*mA,LTraceIndexSetsBR[1],LTraceIndexSetsBR[0]);
  mA01_TR = slice(*mA,LTraceIndexSetsTR[0],LTraceIndexSetsTR[1]);
  mA10_TR = slice(*mA,LTraceIndexSetsTR[1],LTraceIndexSetsTR[0]);
  mA01_TL = slice(*mA,LTraceIndexSetsTL[0],LTraceIndexSetsTL[1]);
  mA10_TL = slice(*mA,LTraceIndexSetsTL[1],LTraceIndexSetsTL[0]);

  mA01Shifted_B = slice(*mA,indexSets2Shifted[1],indexSets2Shifted[2]); 
  mA10Shifted_B = slice(*mA,indexSets2Shifted[2],indexSets2Shifted[1]); 
  mA01Shifted_R = slice(*mA,indexSets1Shifted[5],indexSets1Shifted[4]); 
  mA10Shifted_R = slice(*mA,indexSets1Shifted[4],indexSets1Shifted[5]); 
  mA01Shifted_T = slice(*mA,indexSets2Shifted[5],indexSets2Shifted[4]); 
  mA10Shifted_T = slice(*mA,indexSets2Shifted[4],indexSets2Shifted[5]); 
  mA01Shifted_L = slice(*mA,indexSets1Shifted[1],indexSets1Shifted[2]); 
  mA10Shifted_L = slice(*mA,indexSets1Shifted[2],indexSets1Shifted[1]); 
  mA01Shifted_BL = slice(*mA,LTraceIndexSetsBLShifted[0],LTraceIndexSetsBLShifted[1]);
  mA10Shifted_BL = slice(*mA,LTraceIndexSetsBLShifted[1],LTraceIndexSetsBLShifted[0]);
  mA01Shifted_BR = slice(*mA,LTraceIndexSetsBRShifted[0],LTraceIndexSetsBRShifted[1]);
  mA10Shifted_BR = slice(*mA,LTraceIndexSetsBRShifted[1],LTraceIndexSetsBRShifted[0]);
  mA01Shifted_TR = slice(*mA,LTraceIndexSetsTRShifted[0],LTraceIndexSetsTRShifted[1]);
  mA10Shifted_TR = slice(*mA,LTraceIndexSetsTRShifted[1],LTraceIndexSetsTRShifted[0]);
  mA01Shifted_TL = slice(*mA,LTraceIndexSetsTLShifted[0],LTraceIndexSetsTLShifted[1]);
  mA10Shifted_TL = slice(*mA,LTraceIndexSetsTLShifted[1],LTraceIndexSetsTLShifted[0]);

  mA01ShiftedLR_BL = slice(*mA,LTraceIndexSetsBLShiftedLR[0],LTraceIndexSetsBLShiftedLR[1]);
  mA10ShiftedLR_BL = slice(*mA,LTraceIndexSetsBLShiftedLR[1],LTraceIndexSetsBLShiftedLR[0]);
  mA01ShiftedLR_BR = slice(*mA,LTraceIndexSetsBRShiftedLR[0],LTraceIndexSetsBRShiftedLR[1]);
  mA10ShiftedLR_BR = slice(*mA,LTraceIndexSetsBRShiftedLR[1],LTraceIndexSetsBRShiftedLR[0]);
  mA01ShiftedLR_TR = slice(*mA,LTraceIndexSetsTRShiftedLR[0],LTraceIndexSetsTRShiftedLR[1]);
  mA10ShiftedLR_TR = slice(*mA,LTraceIndexSetsTRShiftedLR[1],LTraceIndexSetsTRShiftedLR[0]);
  mA01ShiftedLR_TL = slice(*mA,LTraceIndexSetsTLShiftedLR[0],LTraceIndexSetsTLShiftedLR[1]);
  mA10ShiftedLR_TL = slice(*mA,LTraceIndexSetsTLShiftedLR[1],LTraceIndexSetsTLShiftedLR[0]);

  mA01ShiftedUD_BL = slice(*mA,LTraceIndexSetsBLShiftedUD[0],LTraceIndexSetsBLShiftedUD[1]);
  mA10ShiftedUD_BL = slice(*mA,LTraceIndexSetsBLShiftedUD[1],LTraceIndexSetsBLShiftedUD[0]);
  mA01ShiftedUD_BR = slice(*mA,LTraceIndexSetsBRShiftedUD[0],LTraceIndexSetsBRShiftedUD[1]);
  mA10ShiftedUD_BR = slice(*mA,LTraceIndexSetsBRShiftedUD[1],LTraceIndexSetsBRShiftedUD[0]);
  mA01ShiftedUD_TR = slice(*mA,LTraceIndexSetsTRShiftedUD[0],LTraceIndexSetsTRShiftedUD[1]);
  mA10ShiftedUD_TR = slice(*mA,LTraceIndexSetsTRShiftedUD[1],LTraceIndexSetsTRShiftedUD[0]);
  mA01ShiftedUD_TL = slice(*mA,LTraceIndexSetsTLShiftedUD[0],LTraceIndexSetsTLShiftedUD[1]);
  mA10ShiftedUD_TL = slice(*mA,LTraceIndexSetsTLShiftedUD[1],LTraceIndexSetsTLShiftedUD[0]);


  // if there is a solver set, set the matrix in the local solver
  if( mSolver != NULL )
  {
    mSolver->setMatrix(mA);
  }
}

// assembles the system matrix,stores it in cahce, and returns it
void SubdomainStructured::assembleSystemMatrix(ComplexSparseMatrix &res ) const
{
  // assemble the system matrix and cache it
  assembleSystemMatrix();

  // copy the system matrix into the result
  res = *mA;
}

// factorizes the liner system
void SubdomainStructured::factorize()
{
  mSolver->factorize();
}

// solves the linear system for an rhs and returns the result
void SubdomainStructured::solve(ComplexVector &x, ComplexVector &res ) const
{
  mSolver->solve(x,res);
}

// returns the mesh-size
Double SubdomainStructured::getMeshSize() const
{
  // if no mesh-size is cached yet, compute it
  if( mH==0.0 )
    mH = mAssembler->getMeshSize();

  // return the mesh size
  return mH;
}

// returns the number of total elements
UInt SubdomainStructured::getNEls() const
{
  return mAssembler->getNEls();
}

// returns the non-shifted or shifted index sets in a parametric direction
void SubdomainStructured::getIndexSets( string traceType, bool shifted, vector<UIntVector> &res ) const
{
  // find the index where the parametric direction is stored
  // if it is not found, return an error that no limits have been set
  auto it = find(mTraceTypes.begin(), mTraceTypes.end(), traceType );
  if( it == mTraceTypes.end() )
  {
    cerr << "ERROR: No limits sets for this parametric direction!" << endl;
    exit(1);
  }
  UInt index = it-mTraceTypes.begin();

  // assign the index sets to the result 
  if( shifted )
    res = mVolumeIndexSetsShifted[index];
  else
    res = mVolumeIndexSets[index];
}

// returns the non-shifted or shifted element index sets in a parametric direction
void SubdomainStructured::getElIndexSets( string traceType, bool shifted, vector<UIntVector> &res ) const
{
  // find the index where the parametric direction is stored
  // if it is not found, return an error that no limits have been set
  auto it = find(mTraceTypes.begin(), mTraceTypes.end(), traceType );
  if( it == mTraceTypes.end() )
  {
    cerr << "ERROR: No limits sets for this parametric direction!" << endl;
    exit(1);
  }
  UInt index = it-mTraceTypes.begin();

  // assign the element index sets to the result 
  if( shifted )
    res = mVolumeElIndexSetsShifted[index];
  else
    res = mVolumeElIndexSets[index];
}

// applies the local system matrix
void SubdomainStructured::applySystemMatrix( const ComplexVector &x, ComplexVector &res ) const
{
  // if the matrix is not assembled yet, assemble it
  if( mA==NULL )
  {
    assembleSystemMatrix();
  }

  // return the result
  res = (*mA)*x;
}

// computes the non-shifted and shifted index sets for a parametric direction
vector<vector<UIntVector>> SubdomainStructured::computeVolumeIndexSets( string traceType, const DoubleVector &lims )
{
  // in this subdomain there are always two parametric direction, since we only support 
  // 2D domain decompositions so far
  vector<vector<UIntVector>> res(2);

  // get the index sets where both the shifted and non-shifted traces are used, this returns a
  // vector with the following index Sets:
  // 0: DOFs below trace0 non-shifted
  // 1: trace0 non-shifted
  // 2: trace1 non-shifted 
  // 3: trace0 shifted
  // 4: trace1 shifted
  // 5: remaining volume dofs
  // 6: traceN non-shifted 
  // 7: traceNP non-shifted 
  // 8: traceN shifted 
  // 9: traceNP shifted 
  //10: DOFS above traceNP shifted
  vector<UIntVector> indexSets = assignIndexSets(traceType,lims);

  // assign the non-shifted index sets
  res[0] = vector<UIntVector>(7);
  res[0][0] = indexSets[0];
  res[0][1] = indexSets[1];
  res[0][2] = indexSets[2];
  res[0][3] = UIntVector(indexSets[3].rows()+indexSets[4].rows()+indexSets[5].rows());
  res[0][3] << indexSets[3],indexSets[4],indexSets[5];
  res[0][4] = indexSets[6];
  res[0][5] = indexSets[7];
  res[0][6] = UIntVector(indexSets[8].rows()+indexSets[9].rows()+indexSets[10].rows());
  res[0][6] << indexSets[8],indexSets[9],indexSets[10];

  // assign the shifted index sets
  res[1] = vector<UIntVector>(7);
  res[1][0] = UIntVector(indexSets[0].rows()+indexSets[1].rows()+indexSets[2].rows());
  res[1][0] << indexSets[0],indexSets[1],indexSets[2];
  res[1][1] = indexSets[3];
  res[1][2] = indexSets[4];
  res[1][3] = UIntVector(indexSets[5].rows()+indexSets[6].rows()+indexSets[7].rows());
  res[1][3] << indexSets[5],indexSets[6],indexSets[7];
  res[1][4] = indexSets[8];
  res[1][5] = indexSets[9];
  res[1][6] = indexSets[10];

  // sort the index sets
  for( UInt i=0; i<res.size(); i++ )
    for( UInt j=0; j<res[i].size(); j++ )
      sort(res[i][j].data(),res[i][j].data()+res[i][j].rows());

  // return the result
  return res;
}

// computes the non-shifted and shifted element index sets for a parametric direction
vector<UIntVector> SubdomainStructured::computeVolumeElIndexSets( string traceType, const vector<UIntVector> &volumeInds ) const
{
  // get the parameter component from the trace type
  UInt component;
  if( traceType=="LeftRight" )
    component = 0;
  else if( traceType=="UpDown" )
    component = 1;
  else
  {
    cerr << "ERROR: Unknown trace type!" << endl;
    exit(1);
  }

  // find the minimum and maximum value of the physcial points in the parametric direction
  vector<DoubleVector> X;
  getPhysPts(X);
  Double minBound = X[0][component];
  Double maxBound = X[0][component];
  for( UInt i=0; i<X.size(); i++ )
  {
    if( X[i][component] > maxBound )
      maxBound = X[i][component];
    if( X[i][component] < minBound )
      minBound = X[i][component];
  }

  // find the minimum parametric value for the four traces
  vector<Double> minLims(4);
  // trace0
  if( volumeInds[1].size()>0 )
    minLims[0]=X[volumeInds[1][0]][component];
  else
    minLims[0]=minBound-1.0;
  for( UInt i=0; i<volumeInds[1].size(); i++ )
    if( X[volumeInds[1][i]][component]<minLims[0] )
      minLims[0] = X[volumeInds[1][i]][component];
  // trace1
  if( volumeInds[2].size()>0 )
    minLims[1]=X[volumeInds[2][0]][component];
  else
    minLims[1]=minBound-1.0;
  for( UInt i=0; i<volumeInds[2].size(); i++ )
    if( X[volumeInds[2][i]][component]<minLims[1] )
      minLims[1] = X[volumeInds[2][i]][component];
  // traceN
  if( volumeInds[4].size()>0 )
    minLims[2]=X[volumeInds[4][0]][component];
  else
    minLims[2]=maxBound+1.0;
  for( UInt i=0; i<volumeInds[4].size(); i++ )
    if( X[volumeInds[4][i]][component]<minLims[2] )
      minLims[2] = X[volumeInds[4][i]][component];
  // traceNP
  if( volumeInds[5].size()>0 )
    minLims[3]=X[volumeInds[5][0]][component];
  else
    minLims[3]=maxBound+1.0;
  for( UInt i=0; i<volumeInds[5].size(); i++ )
    if( X[volumeInds[5][i]][component]<minLims[3] )
      minLims[3] = X[volumeInds[5][i]][component];

  // find the maximum parametric value for the four traces
  vector<Double> maxLims(4);
  // trace0
  if( volumeInds[1].size()>0 )
    maxLims[0]=X[volumeInds[1][0]][component];
  else
    maxLims[0]=minBound-1.0;
  for( UInt i=0; i<volumeInds[1].size(); i++ )
    if( X[volumeInds[1][i]][component]>maxLims[0] )
      maxLims[0] = X[volumeInds[1][i]][component];
  // trace1
  if( volumeInds[2].size()>0 )
    maxLims[1]=X[volumeInds[2][0]][component];
  else
    maxLims[1]=minBound-1.0;
  for( UInt i=0; i<volumeInds[2].size(); i++ )
    if( X[volumeInds[2][i]][component]>maxLims[1] )
      maxLims[1] = X[volumeInds[2][i]][component];
  // traceN
  if( volumeInds[4].size()>0 )
    maxLims[2]=X[volumeInds[4][0]][component];
  else
    maxLims[2]=maxBound+1.0;
  for( UInt i=0; i<volumeInds[4].size(); i++ )
    if( X[volumeInds[4][i]][component]>maxLims[2] )
      maxLims[2] = X[volumeInds[4][i]][component];
  // traceNP
  if( volumeInds[5].size()>0 )
    maxLims[3]=X[volumeInds[5][0]][component];
  else
    maxLims[3]=maxBound+1.0;
  for( UInt i=0; i<volumeInds[5].size(); i++ )
    if( X[volumeInds[5][i]][component]>maxLims[3] )
      maxLims[3] = X[volumeInds[5][i]][component];

  // compute the midpoints of each element
  UInt nEls = getNEls();
  vector<Double> midpoints_XVals(nEls);
  for( UInt i=0; i<nEls; i++ )
  {
    UIntVector elPts;
    getElPhysPtIndices(i, elPts);
    midpoints_XVals[i] = X[elPts[0]][component];
    for( UInt j=1; j<elPts.size(); j++ )
    {
      midpoints_XVals[i] += X[elPts[j]][component];
    }

    midpoints_XVals[i] /= elPts.size();
  }

  // depending on the midpoint locaiton with respect to the
  // minimum/maximum values of the traces, we can determine
  // what element set an element has to belong to
  // we set it accordingly
  // for later use, we also keep track of counts of elements
  // in each set
  vector<UInt> indexSetsInds(nEls);
  vector<UInt> numIndexSetsInds(9);
  for( UInt i=0; i<9; i++ )
    numIndexSetsInds[i] = 0;
  for( UInt i=0; i<nEls; i++ )
  {
    if( midpoints_XVals[i] < minLims[0] )
    {
      indexSetsInds[i] = 0;
      numIndexSetsInds[0]++;
    }
    else if( midpoints_XVals[i] < maxLims[0] )
    {
      indexSetsInds[i] = 1;
      numIndexSetsInds[1]++;
    }
    else if( midpoints_XVals[i] < minLims[1] )
    {
      indexSetsInds[i] = 2;
      numIndexSetsInds[2]++;
    }
    else if( midpoints_XVals[i] < maxLims[1] )
    {
      indexSetsInds[i] = 3;
      numIndexSetsInds[3]++;
    }
    else if( midpoints_XVals[i] < minLims[2] )
    {
      indexSetsInds[i] = 4;
      numIndexSetsInds[4]++;
    }
    else if( midpoints_XVals[i] < maxLims[2] )
    {
      indexSetsInds[i] = 5;
      numIndexSetsInds[5]++;
    }
    else if( midpoints_XVals[i] < minLims[3] )
    {
      indexSetsInds[i] = 6;
      numIndexSetsInds[6]++;
    }
    else if( midpoints_XVals[i] < maxLims[3] )
    {
      indexSetsInds[i] = 7;
      numIndexSetsInds[7]++;
    }
    else
    {
      indexSetsInds[i] = 8;
      numIndexSetsInds[8]++;
    }
  }

  // if there are no DOFS in any of the traces, we just assign all elements to the interior
  if( volumeInds[1].size()==0 && volumeInds[2].size()==0 && volumeInds[4].size()==0 && volumeInds[5].size()==0 )
  {
    for( UInt i=0; i<9; i++ )
      numIndexSetsInds[i]=0;
    numIndexSetsInds[4]=indexSetsInds.size();
    for( UInt i=0; i<indexSetsInds.size(); i++ )
      indexSetsInds[i] = 4;
  }

  // Otherwise, we can construct the element sets
  vector<UIntVector> res(9);
  for( UInt i=0; i<9; i++ )
  {
    res[i] = UIntVector(numIndexSetsInds[i]);
    UInt curr_ind = 0;
    for( UInt j=0; j<indexSetsInds.size(); j++ )
      if( indexSetsInds[j]==i )
      {
        res[i][curr_ind] = j;
        curr_ind++;
      }
  }

  // return the result
  return res;
}

// computes the auxilliary index sets that take the shifted and non-shifted traces into account
// returns the following index sets:
// get the index sets where both the shifted and non-shifted traces are used, this returns a
// vector with the following index Sets:
// 0: DOFs below trace0 non-shifted
// 1: trace0 non-shifted
// 2: trace1 non-shifted 
// 3: trace0 shifted
// 4: trace1 shifted
// 5: remaining volume dofs
// 6: traceN non-shifted 
// 7: traceNP non-shifted 
// 8: traceN shifted 
// 9: traceNP shifted 
//10: DOFS above traceNP shifted
vector<UIntVector> SubdomainStructured::assignIndexSets(string traceType, const DoubleVector &lims)
{
  // lims has to have size 10
  // the limits have to be described in a 10-element vector:
  // Each limit should be in between two rows of degrees of freedom.
  // The limiters are defined to be below (left of) the first row of DOFs
  // corresponding to
  // 0: trace0 non-shifted
  // 1: trace1 non-shifted
  // 2: trace0 shifted
  // 3: trace1 shifted
  // 4: remaining volume DOFs
  // 5: traceN non-shifted
  // 6: traceNP non-shifted
  // 7: traceN shifted
  // 8: traceNP shifted
  // 9: starting row of DOFs above traceNP shifted
  assert( lims.size()==10 );

  // find the parametric direction, if it is not found, print an error and stop 
  auto it = find(mTraceTypes.begin(), mTraceTypes.end(), traceType );
  if( it == mTraceTypes.end() )
  {
    cerr << "ERROR: The trace type before index sets can be computed!" << endl;
    exit(1);
  }
  UInt component = it - mTraceTypes.begin();

  // get the physical points corresponding to the DOFs, and the mesh-size
  vector<DoubleVector> X;
  getDOFPts(X);
  Double h = getMeshSize();

  // assign a domain index to each physical DOF point, the indices have the following meaning:
  // 0: DOFs below trace0 non-shifted
  // 1: trace0 non-shifted
  // 2: trace1 non-shifted 
  // 3: trace0 shifted
  // 4: trace1 shifted
  // 5: remaining volume dofs
  // 6: traceN non-shifted 
  // 7: traceNP non-shifted 
  // 8: traceN shifted 
  // 9: traceNP shifted 
  //10: DOFS above traceNP shifted
  UIntVector domainIndices(X.size());
  for(UInt i=0; i<X.size(); i++ )
  {
    bool done=false;

    for( UInt j=0; j<10; j++ )
    {
      if( X[i][component] < lims[j] )
      {
        domainIndices[i] = j;
        done = true;
        break;
      }

      if( done==false )
        domainIndices[i] = 10;
    }
  }

  // allocate memory for the result
  vector<UIntVector> res(11);
  vector<UInt> numIndex(11);
  vector<UInt> currIndex(11);
  for( UInt i=0; i<11; i++ )
  {
    numIndex[i] = 0;
    currIndex[i] = 0;
  }

  // count how many indices there are in each set
  for( UInt i=0; i<domainIndices.size(); i++ )
  {
    numIndex[domainIndices[i]]++;
  }

  // allocate memory for the result
  for( UInt i=0; i<11; i++ )
    res[i] = UIntVector(numIndex[i]);

  // assign the indices of the index sets
  for( UInt i=0; i<domainIndices.size(); i++ )
  {
    res[domainIndices[i]][currIndex[domainIndices[i]]]=i;
    currIndex[domainIndices[i]]++;
  }

  // sort the index sets
  for( UInt i=0; i<res.size(); i++ )
    sort( res[i].data(), res[i].data()+res[i].size() );

  // return the result
  return res;
}

// the functiont hat computes the L trace indices
vector<UIntVector> SubdomainStructured::computeLTraceIndexSets( string LTraceType, bool shiftedLR, bool shiftedUD )
{
  // get the LeftRight and Updown indices
  vector<UIntVector> indexSets1, indexSets2;
  getIndexSets("LeftRight",shiftedLR,indexSets1);
  getIndexSets("UpDown",shiftedUD,indexSets2);

  // define variables for the vertical and horizontal shifted and non-shfited index sets
  UIntVector ind0_v, ind0_h, ind1_v, ind1_h;
  // define variables for the vertical and horizontal shifted and non-shfited indices that are not part of the L traces
  UIntVector nonInd0_v, nonInd0_h, nonInd1_v, nonInd1_h;

  // for each trace set up the indices
  if( LTraceType=="BottomLeft" )
  {
    // the vertical and horizontal traces are the ones on the left and on the bottom 
    ind0_v = indexSets1[1];    
    ind1_v = indexSets1[2];    
    ind0_h = indexSets2[1];    
    ind1_h = indexSets2[2];   

    // for trace0, the indices on the left are not part of the L trace
    nonInd0_v = indexSets1[0];
    // for trace0, the indices on the bottom are not part of the L trace
    nonInd0_h = indexSets2[0];
    // for trace1, the indices on the left of trace1 in the vertical direction are not part of the L trace
    nonInd1_v = UIntVector(indexSets1[0].size()+indexSets1[1].size());
    nonInd1_v << indexSets1[0],indexSets1[1];
    // for trace1, the indices on the bottom of trace1 in the horizontal direction are not part of the L trace
    nonInd1_h = UIntVector(indexSets2[0].size()+indexSets2[1].size());
    nonInd1_h << indexSets2[0],indexSets2[1];
  }
  else if( LTraceType=="BottomRight" )
  {
    // the vertical and horizontal traces are the ones on the right and on the bottom 
    ind0_v = indexSets1[5];    
    ind1_v = indexSets1[4];    
    ind0_h = indexSets2[1];    
    ind1_h = indexSets2[2];   

    // for trace0, the indices on the right are not part of the L trace
    nonInd0_v = indexSets1[6];
    // for trace0, the indices on the bottom are not part of the L trace
    nonInd0_h = indexSets2[0];
    // for trace1, the indices on the right of trace1 in the vertical direction are not part of the L trace
    nonInd1_v = UIntVector(indexSets1[5].size()+indexSets1[6].size());
    nonInd1_v << indexSets1[5],indexSets1[6];
    // for trace1, the indices on the bottom of trace1 in the horizontal direction are not part of the L trace
    nonInd1_h = UIntVector(indexSets2[0].size()+indexSets2[1].size());
    nonInd1_h << indexSets2[0],indexSets2[1];
  }
  else if( LTraceType=="TopRight" )
  {
    // the vertical and horizontal traces are the ones on the right and on the top
    ind0_v = indexSets1[5];    
    ind1_v = indexSets1[4];    
    ind0_h = indexSets2[5];    
    ind1_h = indexSets2[4];   

    // for trace0, the indices on the right are not part of the L trace
    nonInd0_v = indexSets1[6];
    // for trace0, the indices on the top are not part of the L trace
    nonInd0_h = indexSets2[6];
    // for trace1, the indices on the right of trace1 in the vertical direction are not part of the L trace
    nonInd1_v = UIntVector(indexSets1[5].size()+indexSets1[6].size());
    nonInd1_v << indexSets1[5],indexSets1[6];
    // for trace1, the indices on the top of trace1 in the horizontal direction are not part of the L trace
    nonInd1_h = UIntVector(indexSets2[5].size()+indexSets2[6].size());
    nonInd1_h << indexSets2[5],indexSets2[6];
  }
  else if( LTraceType=="TopLeft" )
  {
    // the vertical and horizontal traces are the ones on the left and on the top
    ind0_v = indexSets1[1];    
    ind1_v = indexSets1[2];    
    ind0_h = indexSets2[5];    
    ind1_h = indexSets2[4];   

    // for trace0, the indices on the left are not part of the L trace
    nonInd0_v = indexSets1[0];
    // for trace0, the indices on the top are not part of the L trace
    nonInd0_h = indexSets2[6];
    // for trace1, the indices on the left of trace1 in the vertical direction are not part of the L trace
    nonInd1_v = UIntVector(indexSets1[0].size()+indexSets1[1].size());
    nonInd1_v << indexSets1[0],indexSets1[1];
    // for trace1, the indices on the top of trace1 in the horizontal direction are not part of the L trace
    nonInd1_h = UIntVector(indexSets2[5].size()+indexSets2[6].size());
    nonInd1_h << indexSets2[5],indexSets2[6];
  }
  else
  {
    cerr << "ERROR: Invalid LTrace type!" << endl;
    exit(1);
  }

  // allocate memory for the result
  vector<UIntVector> res(2);

  // define a temporary variable to set which indices are used to define traces
  vector<bool> temp(getSize());
  for(UInt k=0; k<temp.size();k++ )
    temp[k] = false;

  // set which indices are part of trace0 and not, use the horizontal and vertical traces
  for( UInt k=0; k<ind0_v.size(); k++ )
    temp[ind0_v[k]] = true;
  for( UInt k=0; k<ind0_h.size(); k++ )
    temp[ind0_h[k]] = true;
  for( UInt k=0; k<nonInd0_v.size(); k++ )
    temp[nonInd0_v[k]] = false;
  for( UInt k=0; k<nonInd0_h.size(); k++ )
    temp[nonInd0_h[k]] = false;

  // count how many indices there are in the L trace
  UInt size=0;
  for( UInt k=0; k<temp.size(); k++ )
    if( temp[k] )
      size++;
  // allocate memory and set the trace0 indices
  res[0] = UIntVector(size);
  UInt curr_index=0;
  for( UInt k=0; k<temp.size(); k++ )
    if( temp[k] )
    {
      res[0][curr_index] = k;
      curr_index++;
    }

  // reinitialize temp
  for(UInt k=0; k<temp.size();k++ )
    temp[k] = false;

  // set which indices are part of trace1 and not, use the horizontal and vertical traces
  for( UInt k=0; k<ind1_v.size(); k++ )
    temp[ind1_v[k]] = true;
  for( UInt k=0; k<ind1_h.size(); k++ )
    temp[ind1_h[k]] = true;
  for( UInt k=0; k<nonInd1_v.size(); k++ )
    temp[nonInd1_v[k]] = false;
  for( UInt k=0; k<nonInd1_h.size(); k++ )
    temp[nonInd1_h[k]] = false;

  // count how many indices there are in the L trace
  size=0;
  for( UInt k=0; k<temp.size(); k++ )
    if( temp[k] )
      size++;
  // allocate memory and set the trace0 indices
  res[1] = UIntVector(size);
  curr_index=0;
  for( UInt k=0; k<temp.size(); k++ )
    if( temp[k] )
    {
      res[1][curr_index] = k;
      curr_index++;
    }

  return res;
}

// the function that returns the L trce indices
void SubdomainStructured::getLTraceIndexSets( string LtraceType, bool shiftedLR, bool shiftedUD, vector<UIntVector> &res ) const
{
  // find the index where the LTraceType is stored
  // if it is not found, return an error that no limits have been set
  auto it = find(mLTraceTypes.begin(), mLTraceTypes.end(), LtraceType );
  if( it == mLTraceTypes.end() )
  {
    cerr << "ERROR: No limits sets for this LTraceType!" << endl;
    exit(1);
  }
  UInt index = it-mLTraceTypes.begin();

  // assign the index sets to the result 
  if( shiftedLR && shiftedUD )
    res = mLTraceIndexSetsShifted[index];
  else if( shiftedLR )
    res = mLTraceIndexSetsShiftedLR[index];
  else if( shiftedUD )
    res = mLTraceIndexSetsShiftedUD[index];
  else
    res = mLTraceIndexSets[index];
}

// the function that returns the equivalent sources for a given trace values
void SubdomainStructured::getEquivalentSources( string traceType, bool shiftedLR, bool shiftedUD, ComplexVector vals0, ComplexVector vals1, vector<ComplexVector> &res ) const
{
  // allocate memory for the result
  res = vector<ComplexVector>(2);
  // define pointers to the matrices used for the equivalent sources
  ComplexSparseMatrix *A01, *A10;
  // depending on the traceType define the pointers to the matrices
  if( traceType=="Bottom" )
  {
    if( shiftedUD )
    {
      A01 = &mA01Shifted_B;    
      A10 = &mA10Shifted_B;    
    }
    else
    {
      A01 = &mA01_B;    
      A10 = &mA10_B;    
    }
  }
  else if( traceType=="Right" )
  {
    if( shiftedLR )
    {
      A01 = &mA01Shifted_R;    
      A10 = &mA10Shifted_R;    
    }
    else
    {
      A01 = &mA01_R;    
      A10 = &mA10_R;    
    }
  }
  else if( traceType=="Top" )
  {
    if( shiftedUD )
    {
      A01 = &mA01Shifted_T;    
      A10 = &mA10Shifted_T;    
    }
    else
    {
      A01 = &mA01_T;    
      A10 = &mA10_T;    
    }
  }
  else if( traceType=="Left" )
  {
    if( shiftedLR )
    {
      A01 = &mA01Shifted_L;    
      A10 = &mA10Shifted_L;    
    }
    else
    {
      A01 = &mA01_L;    
      A10 = &mA10_L;    
    }
  }
  else if( traceType=="BottomLeft" )
  {
    if( shiftedLR && shiftedUD )
    {
      A01 = &mA01Shifted_BL;    
      A10 = &mA10Shifted_BL;    
    }
    else if( shiftedLR )
    {
      A01 = &mA01ShiftedLR_BL;    
      A10 = &mA10ShiftedLR_BL;    
    }
    else if( shiftedUD )
    {
      A01 = &mA01ShiftedUD_BL;    
      A10 = &mA10ShiftedUD_BL;    
    }
    else
    {
      A01 = &mA01_BL;    
      A10 = &mA10_BL;    
    }
  }
  else if( traceType=="BottomRight" )
  {
    if( shiftedLR && shiftedUD )
    {
      A01 = &mA01Shifted_BR;    
      A10 = &mA10Shifted_BR;    
    }
    else if( shiftedLR )
    {
      A01 = &mA01ShiftedLR_BR;    
      A10 = &mA10ShiftedLR_BR;    
    }
    else if( shiftedUD )
    {
      A01 = &mA01ShiftedUD_BR;    
      A10 = &mA10ShiftedUD_BR;    
    }
    else
    {
      A01 = &mA01_BR;    
      A10 = &mA10_BR;    
    }
  }
  else if( traceType=="TopRight" )
  {
    if( shiftedLR && shiftedUD )
    {
      A01 = &mA01Shifted_TR;    
      A10 = &mA10Shifted_TR;    
    }
    else if( shiftedLR )
    {
      A01 = &mA01ShiftedLR_TR;    
      A10 = &mA10ShiftedLR_TR;    
    }
    else if( shiftedUD )
    {
      A01 = &mA01ShiftedUD_TR;    
      A10 = &mA10ShiftedUD_TR;    
    }
    else
    {
      A01 = &mA01_TR;    
      A10 = &mA10_TR;    
    }
  }
  else if( traceType=="TopLeft" )
  {
    if( shiftedLR && shiftedUD )
    {
      A01 = &mA01Shifted_TL;    
      A10 = &mA10Shifted_TL;    
    }
    else if( shiftedLR )
    {
      A01 = &mA01ShiftedLR_TL;    
      A10 = &mA10ShiftedLR_TL;    
    }
    else if( shiftedUD )
    {
      A01 = &mA01ShiftedUD_TL;    
      A10 = &mA10ShiftedUD_TL;    
    }
    else
    {
      A01 = &mA01_TL;    
      A10 = &mA10_TL;    
    }
  }

  // check consistencies
  assert( vals0.size()==A01->rows() );
  assert( vals1.size()==A01->cols() );
  assert( vals1.size()==A10->rows() );
  assert( vals0.size()==A10->cols() );

  // define equivalent sources
  res[0] = (*A01)*vals1;
  res[1] = -(*A10)*vals0;
}

// the function that extracts L trace values fro horizontal and vertiacal trace values
void SubdomainStructured::extractLTraces( string LTraceType, bool shiftedLR, bool shiftedUD, ComplexVector valsLR0, ComplexVector valsLR1, ComplexVector valsUD0, ComplexVector valsUD1, vector<ComplexVector> &res ) const
{
  // define a vector of values that is as big as the subdomain size
  ComplexVector tempVals(getSize());

  // get the horizontal and vertical trace indices
  vector<UIntVector> indexSets1, indexSets2;
  getIndexSets("LeftRight",shiftedLR, indexSets1);
  getIndexSets("UpDown",shiftedUD, indexSets2);

  // define variables for the horizontal and vertical trace indices
  UIntVector indLR0, indLR1, indUD0, indUD1;
  // depending on the L trace type, define what is used for the horizontal and vertialc trces
  if( LTraceType=="BottomLeft" )
  {
    indLR0 = indexSets1[1]; 
    indLR1 = indexSets1[2]; 
    indUD0 = indexSets2[1]; 
    indUD1 = indexSets2[2]; 
  }
  else if( LTraceType=="BottomRight" )
  {
    indLR0 = indexSets1[5];
    indLR1 = indexSets1[4];
    indUD0 = indexSets2[1]; 
    indUD1 = indexSets2[2]; 
  }
  else if( LTraceType=="TopRight" )
  {
    indLR0 = indexSets1[5];
    indLR1 = indexSets1[4];
    indUD0 = indexSets2[5]; 
    indUD1 = indexSets2[4]; 
  }
  else if( LTraceType=="TopLeft" )
  {
    indLR0 = indexSets1[1];
    indLR1 = indexSets1[2];
    indUD0 = indexSets2[5]; 
    indUD1 = indexSets2[4]; 
  }

  // use the horizontal and vertical trace indices and values to set the values in the temp variable
  assert(valsLR1.size()==indLR1.size());
  for( UInt k=0; k<indLR1.size(); k++ )
    tempVals[indLR1[k]]=valsLR1[k];

  assert(valsUD1.size()==indUD1.size());
  for( UInt k=0; k<indUD1.size(); k++ )
    tempVals[indUD1[k]]=valsUD1[k];

  assert(valsLR0.size()==indLR0.size());
  for( UInt k=0; k<indLR0.size(); k++ )
    tempVals[indLR0[k]]=valsLR0[k];

  assert(valsUD0.size()==indUD0.size());
  for( UInt k=0; k<indUD0.size(); k++ )
    tempVals[indUD0[k]]=valsUD0[k];

  // use the temp variable to extract the L traces
  vector<UIntVector> traceInds;
  // get the L trace indices
  getLTraceIndexSets( LTraceType, shiftedLR, shiftedUD, traceInds );
  // allocate memory and set the values using the L trace indices and the temp variable
  res = vector<ComplexVector>(2);
  res[0] = ComplexVector(traceInds[0].size());
  for( UInt k=0; k<traceInds[0].size(); k++ )
    res[0][k] = tempVals[traceInds[0][k]];
  res[1] = ComplexVector(traceInds[1].size());
  for( UInt k=0; k<traceInds[1].size(); k++ )
    res[1][k] = tempVals[traceInds[1][k]];
}

