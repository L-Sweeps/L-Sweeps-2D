#include "RHSSrcPts.h"
#include <iostream>
#include "utils.h"

using namespace std;

// The standard constructor, sets up an empty object
RHSSrcPts::RHSSrcPts(): RHS(), mEps(0)
{

}

// The cop constructor
RHSSrcPts::RHSSrcPts( const RHSSrcPts &other ): RHS(other), mSourcePts( other.mSourcePts ), mEps( other.mEps ) 
{

}

// The destructor
RHSSrcPts::~RHSSrcPts()
{

}

// The assignment operator
RHSSrcPts &RHSSrcPts::operator=( const RHSSrcPts &other )
{
  // only copy if other is different from this object
  if( this != &other )
  {
    // copy all member variables
    mSourcePts = other.mSourcePts;
    mEps = other.mEps;
  }

  return *this;
}

// The function that returns a pointer to a copy
RHS *RHSSrcPts::copy() const
{
  return new RHSSrcPts(*this);
}

// The function clearing out the source point locations
void RHSSrcPts::clearSrcPts()
{
  mSourcePts.clear();
}

// The function that sets the accuracy of approximation by exponential functions
void RHSSrcPts::setEpsilon( Double epsilon )
{
  mEps = epsilon;
}

// The function that sets the source point locations
void RHSSrcPts::setSrcPts( const vector<DoubleVector> &sourcePts )
{
  mSourcePts = sourcePts;
}

// The functin that adds a point source location
void RHSSrcPts::addSrcPt( const DoubleVector &src )
{
  if( mSourcePts.size()>0 && mSourcePts[0].size() != src.size() )
  {
    cerr << "ERROR: Only source points with the same dimension can be added" << endl;
    exit(1);
  }
  mSourcePts.push_back( src );
}

// The function that returns the point source locations
void RHSSrcPts::getSrcPts( vector<DoubleVector> &res ) const
{
  res = mSourcePts;
}

// The funciton that returns the current accuracy to which the point sources are approximated
Double RHSSrcPts::getEpsilon() const
{
  return mEps;
}

// The function that evaluates the right hand side at a point x
Complex RHSSrcPts::eval( DoubleVector x ) const
{
  // Allocate memory for the result
  Complex res = 0;
  
  // loop over all point sources, evaluate the approximation by an exponential function and add the value to the result
  for( UInt i=0; i<mSourcePts.size(); i++ )
    res += deltaFunction( mSourcePts[i], x );

  // return the result
  return res;
}

// The function that approximates a point source by an exponential funcion
Double RHSSrcPts::deltaFunction( const DoubleVector &x0, const DoubleVector &x ) const
{
  // check if the input is consistend and get the accuracy
  assert( x.size()==x0.size() );
  Double h = mEps;

  // evaluate the approximation and return the result
  Double absXX0sqrd = 0;
  for( UInt i=0; i<x.size(); i++ )
    absXX0sqrd += pow( x0[i]-x[i],2.0 );
  return 1/(pi*pow(h,x.size()))*exp(-absXX0sqrd/(pow(h,2.0)));
}


