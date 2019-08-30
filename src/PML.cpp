#include "PML.h"
#include <iostream>

using namespace::std;

PML::PML(): mParent(NULL), mCompDomain(0), mIntDomain(0), mPMLWidth(0.0)
{

}

// the implementation of the copy constructor
PML::PML( const PML &other ): mParent(other.mParent), mCompDomain(other.mCompDomain), mIntDomain(other.mIntDomain), mPMLWidth(other.mPMLWidth)
{

}

// the implementation of the assignment operator
PML &PML::operator=( const PML &other )
{
  if( &other != this )
  {
    mParent = other.mParent;
    mCompDomain = other.mCompDomain;
    mIntDomain = other.mIntDomain;
    mPMLWidth = other.mPMLWidth;
  }
  return *this;
}

void PML::setParent( const AssemblerHelmholtz *parent )
{
  mParent = parent;
}

void PML::setWidth( Double pmlWidth )
{
  mPMLWidth = pmlWidth;
  if( mIntDomain.size()>0 )
  {
    mCompDomain = mIntDomain;
    mCompDomain[0][0] -= pmlWidth;
    mCompDomain[0][1] += pmlWidth;
    mCompDomain[1][0] -= pmlWidth;
    mCompDomain[1][1] += pmlWidth;
  }
}

Double PML::getWidth() const
{
  return mPMLWidth;
}

void PML::setDomain( const vector<DoubleVector> &domain )
{
  mIntDomain = domain;
  if( mPMLWidth>0 )
  {
    mCompDomain = mIntDomain;
    mCompDomain[0][0] -= mPMLWidth;
    mCompDomain[0][1] += mPMLWidth;
    mCompDomain[1][0] -= mPMLWidth;
    mCompDomain[1][1] += mPMLWidth;
  }
}

vector<DoubleVector> PML::getDomain() const
{
  return mIntDomain;
}

