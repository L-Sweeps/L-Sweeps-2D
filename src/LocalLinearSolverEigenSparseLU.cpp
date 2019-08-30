#include <Eigen/SparseLU>
#include "LocalLinearSolverEigenSparseLU.h"
#include <iostream>

using namespace::Eigen;
using namespace::std;

// implementation of the (private) copy constructor
LocalLinearSolverEigenSparseLU::LocalLinearSolverEigenSparseLU( const LocalLinearSolverEigenSparseLU &other ): mIsFactorized(other.mIsFactorized)
{
  mA = other.mA;
}
    
// implementation of the (private) copy constructor
LocalLinearSolverEigenSparseLU &LocalLinearSolverEigenSparseLU::operator=(const LocalLinearSolverEigenSparseLU &other)
{
  if( this != &other )
  {
    mA = other.mA;
    mIsFactorized = other.mIsFactorized;
  }
  return *this;
}

// generate a copy of this object and return a point to it
LocalLinearSolver *LocalLinearSolverEigenSparseLU::copy() const
{
  return new LocalLinearSolverEigenSparseLU(*this);
}

// sets the matrix of the linear system, the ComplexSparse matrix
// has to be compressed before this function is called, this is inherited from Eigen
void LocalLinearSolverEigenSparseLU::setMatrix( ComplexSparseMatrix *A )
{
  // set the matrix
  mA = A;
}

// factorizes the system
void LocalLinearSolverEigenSparseLU::factorize()
{
  if( mIsFactorized )
  {
    return;
  }
  // make sure that a matrix is set
  assert( mA!=NULL );

  // factorize using Eigen
  mSolver.compute(*mA);
  
  // check if the factorization was successful
  // if not, print an error and exit
  if( mSolver.info()!=0 )
  {
    cerr << "ERROR: Matrix was not successfully factorized!" << endl;
    exit(1);
  }

  // set the flag that the system is factorized
  mIsFactorized = true;
}

// solves the linear system for the rhs b and returns the result in res
void LocalLinearSolverEigenSparseLU::solve(ComplexVector &b, ComplexVector &res) const
{
  // if the system is not factorized, do it
  if( !mIsFactorized )
  {
    cerr << "ERROR: The linear system has to be factorized before it can be solved!" << endl;
    exit(1);
  }

  // solve the system using Eigen
  res = mSolver.solve(b);
}
