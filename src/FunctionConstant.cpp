#include "FunctionConstant.h"
#include <iostream>

// The standard constructor
FunctionConstant::FunctionConstant()
{

}

// The copz constructor
FunctionConstant::FunctionConstant( const FunctionConstant &other ): mVal(other.mVal), mDom(other.mDom)
{

}

// The destructor
FunctionConstant::~FunctionConstant()
{

}

// The assignment operator
FunctionConstant &FunctionConstant::operator=( const FunctionConstant &other )
{
  // Only copy if other is different from this object
  if( this != &other )
  {
    // copy all member variables
    mVal = other.mVal;
    mDom = other.mDom;
  }

  return *this;
}

//The function clearing the cache, nothing to be done here
void FunctionConstant::clearCache() const
{

}

// The function setting the domain of definition
void FunctionConstant::setDomain( const vector<DoubleVector> &dom ) const
{
    mDom=dom;
}

// The function returning the domain of definition
void FunctionConstant::getDomain( vector<DoubleVector> &res ) const
{
    res = mDom;
}

// The function that sets the value of the constant function
void FunctionConstant::setVal( Complex val )
{
    mVal = val;
}

// The copy function to properly copy the object
Function *FunctionConstant::copy() const
{
    return new FunctionConstant(*this);
}

// The evaluation funtion, only the value has to be returned
Complex FunctionConstant::eval( DoubleVector x ) const
{
    return mVal;
}

// The function returning the minimum real part of the value (i.e. just the real part of the value)
Double FunctionConstant::getMinReal() const
{
  return real(mVal);
}

// The function returning the maximum real part of the value (i.e. just the real part of the value)
Double FunctionConstant::getMaxReal() const
{
  return real(mVal);
}

// The function returning the minimum imaginary part of the value (i.e. just the real part of the value)
Double FunctionConstant::getMinImag() const
{
  return imag(mVal);
}

// The function returning the maximum imaginary part of the value (i.e. just the real part of the value)
Double FunctionConstant::getMaxImag() const
{
  return imag(mVal);
}

