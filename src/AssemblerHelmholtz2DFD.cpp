#include "AssemblerHelmholtz2DFD.h"
#include <iostream>
#include "FiniteDifferences.h"
#include "utils.h"
#include "RHS.h"
#include "Function.h"

using namespace::std;

// sets up an empty object
AssemblerHelmholtz2DFD::AssemblerHelmholtz2DFD(): AssemblerHelmholtz(), mOrder(0), mN(0), mIntDomain(0), mCompDomain(0), mPts(0)
{

}

// implementation of the copy constructor
AssemblerHelmholtz2DFD::AssemblerHelmholtz2DFD( const AssemblerHelmholtz2DFD &other ): AssemblerHelmholtz( other ), mN( other.mN ), mOrder( other.mOrder ), mIntDomain(other.mIntDomain), mCompDomain(other.mCompDomain), mPts(other.mPts)
{

}

// implementation of the assignment operator
AssemblerHelmholtz2DFD &AssemblerHelmholtz2DFD::operator=( const AssemblerHelmholtz2DFD &other )
{
  if( this != &other )
  {
    // copy as AssemblerHelmholtz object
    (AssemblerHelmholtz &)(*this) = (const AssemblerHelmholtz &)other;

    // copy the other objects
    mN = other.mN;
    mOrder = other.mOrder;
    mPts = other.mPts;
    mIntDomain = other.mIntDomain;
    mCompDomain = other.mCompDomain;
  }
  return *this;
}

// generates a proper copy and returns a pointer to it
Assembler *AssemblerHelmholtz2DFD::copy() const
{
  return new AssemblerHelmholtz2DFD(*this);
}

// sets the domain of the assembler
void AssemblerHelmholtz2DFD::setDomain( const vector<DoubleVector> &domain )
{
  // set the interior domain
  mIntDomain = domain;

  // set the domain for the PML
  mPML->setDomain(domain);

  // clear cache and the number of elements in one direction
  mPts = vector<DoubleVector>(0);
  mN = UIntVector(0);
}

// returns the interior domain
vector<DoubleVector> AssemblerHelmholtz2DFD::getDomain() const
{
  return mIntDomain;
}

// sets mN (the number of elements in each direction)
void AssemblerHelmholtz2DFD::setN( const UIntVector &N )
{
  // check if an interior domain is set, it has to be set first
  if( mIntDomain.size()==0 )
  {
    cerr << "ERROR: Domain has to be set before setting N" << endl;
    exit(1);
  }

  // compute the mesh sizes in each direction
  Double h1 = (mIntDomain[0][1]-mIntDomain[0][0])/N[0];
  Double h2 = (mIntDomain[1][1]-mIntDomain[1][0])/N[1];
  // define the global mesh size as the smallest in eahc direction
  Double h = min(h1,h2);

  // adjust the PML width so that it is a multiple of the mesh size
  // get the PML width
  Double pmlWidth = mPML->getWidth();

  // determin the minimum number of elements needed to at least get
  // a widht of pmlWidth using a multiple of the mesh size
  Double pmlWidthDivH = pmlWidth/h;
  UInt pmlEls = pmlWidth/h;
  if( pmlWidthDivH-pmlEls>0.5 )
    pmlEls+=1;

  // compute the new width
  pmlWidth = h*pmlEls;

  // set the PML width and the domain
  mPML->setWidth(pmlWidth);
  mPML->setDomain(mIntDomain);
  
  // use the new pml width to define the computational domain
  mCompDomain = mIntDomain;
  mCompDomain[0][0] -= pmlWidth;
  mCompDomain[0][1] += pmlWidth;
  mCompDomain[1][0] -= pmlWidth;
  mCompDomain[1][1] += pmlWidth;

  // compute the number of elements in each direction using the computational domain
  mN = UIntVector(2);
  mN[0] = (mCompDomain[0][1]-mCompDomain[0][0]+0.5*h1)/h1;
  mN[1] = (mCompDomain[1][1]-mCompDomain[1][0]+0.5*h2)/h2;

  // clear the cache of the DOF points
  mPts = vector<DoubleVector>(0);
}

// returns the number of elements in each direction
UIntVector AssemblerHelmholtz2DFD::getN() const
{
  return mN;
}

// sets the order
void AssemblerHelmholtz2DFD::setOrder( UInt order )
{
  mOrder = order;
}

// return the order
UInt AssemblerHelmholtz2DFD::getOrder() const
{
  return mOrder;
}

// returns the mesh-size
Double AssemblerHelmholtz2DFD::getMeshSize() const
{
  // compute the mesh size in the x-direction
  Double h1 = (mCompDomain[0][1]-mCompDomain[0][0])/mN[0];
  // compute the mesh size in the y-direction
  Double h2 = (mCompDomain[1][1]-mCompDomain[1][0])/mN[1];

  // return the mesh-size
  return sqrt(h1*h2);
}

// returns the physical points associated with the DOFs
void AssemblerHelmholtz2DFD::getDOFPts( vector<DoubleVector> &res ) const
{
  // compute the number of DOF points if they are not set yet
  if( mPts.size()==0 )
  {
    // get the number of elements in each direction
    UInt nx = mN[0];
    UInt ny = mN[1];

    // compute the mesh-sizes in each direction
    Double hx = (mCompDomain[0][1]-mCompDomain[0][0])/nx;
    Double hy = (mCompDomain[1][1]-mCompDomain[1][0])/ny;

    // allocate memory for the result
    // there are (nx-1)*(ny-1) physical points because the boundary points
    // are not counted as DOFs, they are just set to zero (due to the PMLs)
    mPts = vector<DoubleVector>((nx-1)*(ny-1));

    // loop over all DOFS and compute their physical location
    for( UInt i=0; i<nx-1; i++)
    {
      for( UInt j=0; j<ny-1; j++)
      {
        // the index of the DOF
        UInt index = i*(ny-1)+j;

        // allocate memory and assign the locations
        mPts[index] = DoubleVector(2);
        mPts[index][0] = mCompDomain[0][0]+(i+1)*hx;
        mPts[index][1] = mCompDomain[1][0]+(j+1)*hy;
      }
    }
  }
  res = mPts;
}

// get the physical ponts that are not DOFs
// (the boundary ponts excluded to set the boundary conditions)
void AssemblerHelmholtz2DFD::getNonDOFPhysPts( vector<DoubleVector> &res ) const
{
  // get the number of elements in each direction
  UInt nx = mN[0];
  UInt ny = mN[1];

  // compute the mesh-sizes in each direction
  Double hx = (mCompDomain[0][1]-mCompDomain[0][0])/nx;
  Double hy = (mCompDomain[1][1]-mCompDomain[1][0])/ny;

  // allocate memory for the result
  res = vector<DoubleVector>(2*nx+2*ny);
  // compute the points along the bottom boundary
  for( UInt i=0; i<nx; i++ )
  {
    res[i] = DoubleVector(2);
    res[i][0] = mCompDomain[0][0]+i*hx;
    res[i][1] = mCompDomain[1][0];
  }
  // compute the points along the right boundary
  for( UInt i=0; i<ny; i++ )
  {
    res[nx+i] = DoubleVector(2);
    res[nx+i][0] = mCompDomain[0][1];
    res[nx+i][1] = mCompDomain[1][0]+hy*i;
  }
  // compute the points along the top boundary
  for( UInt i=0; i<nx; i++ )
  {
    res[nx+ny+i] = DoubleVector(2);
    res[nx+ny+i][0] = mCompDomain[0][1]-i*hx;
    res[nx+ny+i][1] = mCompDomain[1][1];
  }
  // compute the points along the left boundary
  for( UInt i=0; i<ny; i++ )
  {
    res[2*nx+ny+i] = DoubleVector(2);
    res[2*nx+ny+i][0] = mCompDomain[0][0];
    res[2*nx+ny+i][1] = mCompDomain[1][1]-hy*i;
  }
}

// get the assignment of element to physical ponts
// the points are ordered in counter-clockwise order along the boundary of the element
void AssemblerHelmholtz2DFD::getElPhysPtIndices( UInt index, UIntVector &res ) const
{
  // check if the index is consistent
  assert( index<mN[0]*mN[1] );

  // compute the row- and column-index of the element
  UInt J = index/mN[1];
  UInt I = index-J*mN[1];

  // the total number of DOF points
  // this number is important becuase here we take all physical points into
  // account and the physical ponts are ordered as DOF points first and 
  // nonDOF points last
  UInt size = (mN[0]-1)*(mN[1]-1);

  // check if the element indices are consistent
  assert( I<mN[1] );
  assert( J<mN[0] );
  // assign the indices for the bottom left corner element
  if( J==0 && I==0 )
  {
    res = UIntVector(4);
    res[0] = size;
    res[1] = size+1;
    res[2] = 0;
    res[3] = size+2*mN[0]+2*mN[1]-1;
  }
  // assign the indices for the bottom right corner element
  else if( J==mN[0]-1 && I==0 )
  {
    res = UIntVector(4);
    res[0] = size+mN[0]-1;
    res[1] = size+mN[0];
    res[2] = size+mN[0]+1;
    res[3] = (mN[0]-2)*(mN[1]-1);
  }
  // assign the indices for the top right corner element
  else if( J==mN[0]-1 && I==mN[1]-1 )
  {
    res = UIntVector(4);
    res[0] = (mN[0]-1)*(mN[1]-1)-1;
    res[1] = size+mN[0]+mN[1]-1;
    res[2] = size+mN[0]+mN[1];
    res[3] = size+mN[0]+mN[1]+1;
  }
  // assign the indices for the top left corner element
  else if( J==0 && I==mN[1]-1 )
  {
    res = UIntVector(4);
    res[0] = size+2*mN[0]+mN[1]+1;
    res[1] = mN[1]-2;
    res[2] = size+2*mN[0]+mN[1]-1;
    res[3] = size+2*mN[0]+mN[1];
  }
  // assign the indices for an element on the bottom boundary 
  else if( I==0 )
  {
    res = UIntVector(4);
    res[0] = size+J;
    res[1] = size+J+1;
    res[2] = J*(mN[1]-1);
    res[3] = (J-1)*(mN[1]-1);
  }
  // assign the indices for an element on the right boundary 
  else if( J==mN[0]-1 )
  {
    res = UIntVector(4);
    res[0] = (mN[0]-2)*(mN[1]-1)+I-1;
    res[1] = size+mN[0]+I; 
    res[2] = size+mN[0]+I+1; 
    res[3] = (mN[0]-2)*(mN[1]-1)+I;
  }
  // assign the indices for an element on the top boundary 
  else if( I==mN[1]-1 )
  {
    res = UIntVector(4);
    res[0] = (J-1)*(mN[1]-1)+mN[1]-2;
    res[1] = J*(mN[1]-1)+mN[1]-2;
    res[2] = size+2*mN[0]+mN[1]-2-(J-1);
    res[3] = size+2*mN[0]+mN[1]-1-(J-1);
  }
  // assign the indices for an element on the left boundary 
  else if( J==0 )
  {
    res = UIntVector(4);
    res[0] = size+2*mN[0]+2*mN[1]-I;
    res[1] = I-1; 
    res[2] = I;
    res[3] = size+2*mN[0]+2*mN[1]-1-I;
  }
  // assign the indices for an interior element
  else
  {
    res = UIntVector(4);
    res[0] = (J-1)*(mN[1]-1)+(I-1); 
    res[1] = (J-1)*(mN[1]-1)+I; 
    res[2] = J*(mN[1]-1)+I; 
    res[3] = J*(mN[1]-1)+(I-1); 
  }
}

// expands a vector of values corresponding to DOF points to all physical points
// the extension is done by zero
void AssemblerHelmholtz2DFD::expandDOFValsToPhysVals( ComplexVector &res ) const
{
  // allocate memory
  ComplexVector temp((mN[0]+1)*(mN[1]+1));
  // set all values to zero
  for( UInt i=0; i<temp.size(); i++ )
    temp[i] = 0.0;
  // set the values of corresponding to the DOF points
  for( UInt i=0; i<res.size(); i++ )
    temp[i] = res[i];
  // set the new result
  res = temp;
}

// assembles the rhs for a given source function 
void AssemblerHelmholtz2DFD::assembleRHS( const RHS *rhs, ComplexVector &res ) const
{
  // get the physical DOF points
  vector<DoubleVector> pts;
  getDOFPts( pts );

  // allocate memory for the result and initialize to zero
  res = ComplexVector(pts.size());
  for( UInt i=0; i<pts.size(); i++ )
    res[i] = 0;

  // loop over all source points and evaluate the source function at the physical
  // DOF points
  // loop over physical DOF points
  for( UInt i=0; i<pts.size(); i++ )
  {
    // get the PML coefficients for the rhs
    Complex c = mPML->getSourceCoefficient(pts[i]);
    // add the contribution from this source point to the rhs at the DOF pont
    res[i] += c*rhs->eval( pts[i] );
  }
}

// assembles the system matrix
void AssemblerHelmholtz2DFD::assembleSystemMatrix( ComplexSparseMatrix &res ) const
{
  // get the number of elements and the mesh-sizes 
  UInt nx = mN[0];
  UInt ny = mN[1];
  double hx = (mCompDomain[0][1]-mCompDomain[0][0])/nx;
  double hy = (mCompDomain[1][1]-mCompDomain[1][0])/ny;

  // the matrix discretizing the second derivitive in the x-direction
  ComplexSparseMatrix A1 = DerivativeDifferenceMatrix1D( nx+1, hx, mOrder, 2 ).block(1,1,nx-1,nx-1).cast<Complex>();
  
  // the matrix discretizing the second derivitive in the y-direction
  ComplexSparseMatrix A2 = DerivativeDifferenceMatrix1D( ny+1, hy, mOrder, 2 ).block(1,1,ny-1,ny-1).cast<Complex>();
  // the identity matrices in each direction
  ComplexSparseMatrix I1(A1.rows(), A1.cols());
  ComplexSparseMatrix I2(A2.rows(), A2.cols());
  I1.setIdentity();
  I2.setIdentity();

  // get the physical positions of the DOFs
  vector<DoubleVector> pts;
  getDOFPts( pts );

  // compute the kronecker product of A1 and I2
  // this realizes the second derivitive in the x-direction
  A1 = kron(A1,I2);
  assert(pts.size()==A1.rows());
  // compute the kronecker product of I1 and A2
  // this realizes the second derivitive in the y-direction
  A2 = kron(I1,A2);
  assert(pts.size()==A2.rows());

  ComplexVector Lambda1(I2.rows()*I1.rows());
  ComplexVector Lambda2(I2.rows()*I1.rows());
  // loop over the DOF points and assign the values
  for( UInt i=0; i<A1.rows(); i++ )
  {
    ComplexVector Lambda = mPML->getDiffusionCoefficient(pts[i]);
    Lambda1(i) = Lambda[0]; 
    Lambda2(i) = Lambda[1]; 
  }
  // use the PML values and the derivative matrices to realize the term involving
  // the second derivative
  res = -(Lambda1.asDiagonal()*A1)-Lambda2.asDiagonal()*A2;

  // the matrix discretizing the first derivitive in the x-direction
  A1 = DerivativeDifferenceMatrix1D( nx+1, hx, mOrder, 1 ).block(1,1,nx-1,nx-1).cast<Complex>();
  // the matrix discretizing the first derivitive in the y-direction
  A2 = DerivativeDifferenceMatrix1D( ny+1, hy, mOrder, 1 ).block(1,1,ny-1,ny-1).cast<Complex>();
  // compute the kronecker product of A1 and I2
  // this realizes the first derivitive in the x-direction
  A1 = kron(A1,I2);
  assert(pts.size()==A1.rows());
  // compute the kronecker product of I1 and A2
  // this realizes the first derivitive in the y-direction
  A2 = kron(I1,A2);
  assert(pts.size()==A2.rows());
  
  ComplexVector beta1(I2.rows()*I1.rows());
  ComplexVector beta2(I2.rows()*I1.rows());
  // loop over the DOF points and assign the values
  for( UInt i=0; i<A1.rows(); i++ )
  {
    ComplexVector beta = mPML->getConvectionCoefficient(pts[i]);
    beta1(i) = beta[0]; 
    beta2(i) = beta[1]; 
  }
  // use the PML values and the derivative matrices to realize the term involving
  // the second derivative
  res -= beta1.asDiagonal()*A1+beta2.asDiagonal()*A2;

  ComplexVector A1vec(I2.rows()*I1.rows());
  // loop over the DOF points and assign the values
  for( UInt i=0; i<A1.rows(); i++ )
  {
    Complex c = mPML->getReactionCoefficient(pts[i]);
    Complex mVal = mM->eval(pts[i]);
    A1vec(i) = mOmega*mOmega*mVal*c;
  }
  // add to the result
  res -= A1vec.asDiagonal();

  // clear the cache of the wave speed function
  // This means that in subsequent assemblies, the wave speed
  // has to be set again becuase it will be cleard out.
  mM->clearCache();
}

// evaluates an approximation of the delta function centered at x0
// the evaluation point is x
Double AssemblerHelmholtz2DFD::deltaFunction( DoubleVector x0, DoubleVector x ) const
{
  // get the mesh size
  Double h = getMeshSize();
  // compute and return the approximation
   return 1/(pi*h*h)*exp(-(pow(x0[0]-x[0],2.0)+pow(x0[1]-x[1],2.0))/(h*h));
}

// returns the total number of elements
UInt AssemblerHelmholtz2DFD::getNEls() const
{
  return mN[0]*mN[1];
}

