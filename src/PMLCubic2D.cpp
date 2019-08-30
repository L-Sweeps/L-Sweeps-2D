#include <iostream>
#include "PMLCubic2D.h"
#include "utils.h"

// the implementation of the copy constructor
PMLCubic2D::PMLCubic2D( const PMLCubic2D &other ): PML(other), mC(other.mC)
{

}

// the implementation of the copy function
PML *PMLCubic2D::copy() const
{
  return (PML*)new PMLCubic2D(*this);
}

// the implementation of the assignment operator
PMLCubic2D &PMLCubic2D::operator=( const PMLCubic2D &other )
{
  if( &other != this )
  {
    (PML&)*this = (PML&)(other);
    mC = other.mC;
  }
  return *this;
}

// the function to set the absorption coefficient
void PMLCubic2D::setAbsorptionCoeff( Double C )
{
  mC = C;
}

// the function that returns the PML factor for the diffusion term
ComplexVector PMLCubic2D::getDiffusionCoefficient( DoubleVector x ) const
{
  ComplexVector res(2);
  res[0] = pow(alpha1(x),2.0);
  res[1] = pow(alpha2(x),2.0);

  return res;
}

// the function that returns the PML factor for the convection term
ComplexVector PMLCubic2D::getConvectionCoefficient( DoubleVector x ) const
{
  ComplexVector res(2);
  res[0] = alpha1Dx(x)*alpha1(x);
  res[1] = alpha2Dy(x)*alpha2(x);

  return res;
}

// the function that returns the PML factor for the reaction term
Complex PMLCubic2D::getReactionCoefficient( DoubleVector x ) const
{
  return 1.0;
}

// the function that returns the PML factor for the source term
Complex PMLCubic2D::getSourceCoefficient( DoubleVector x ) const
{
  return 1.0;
}

// the auxiliary function alpha1 in the PML 
Complex PMLCubic2D::alpha1( DoubleVector x ) const
{
  // check if the input is good
  assert( x.size()==2 );
  // defime the imaginary unit
  Complex im(0.0,1.0);
  // get omega
  Double omega = mParent->getOmega();
  // return the value
  return 1.0/(1.0+im*sigma1(x,omega));
}

// the derivative of alpha1
Complex PMLCubic2D::alpha1Dx( DoubleVector x ) const
{
  // check if the input is good
  assert( x.size()==2 );
  // defime the imaginary unit
  Complex im(0.0,1.0);
  // get omega
  Double omega = mParent->getOmega();
  // return the value
  return -im*sigma1Dx(x,omega)/(pow(1.0+im*sigma1(x,omega),2.0));
}

// the auxiliary function alpha2 in the PML 
Complex PMLCubic2D::alpha2( DoubleVector x ) const
{
  // check if the input is good
  assert( x.size()==2 );
  // defime the imaginary unit
  Complex im(0.0,1.0);
  // get omega
  Double omega = mParent->getOmega();
  // return the value
  return 1.0/(1.0+im*sigma2(x,omega));
}

// the derivative of alpha2
Complex PMLCubic2D::alpha2Dy( DoubleVector x ) const
{
  // check if the input is good
  assert( x.size()==2 );
  // defime the imaginary unit
  Complex im(0.0,1.0);
  // get omega
  Double omega = mParent->getOmega();
  // return the value
  return -im*sigma2Dy(x,omega)/(pow(1.0+im*sigma2(x,omega),2.0));
}

// the auxiliary function sigma1 in the PML 
Double PMLCubic2D::sigma1( DoubleVector x, Double omega ) const
{
  // check if the input is good
  assert( x.size()==2 );

  // optain the x value
  Double xVal = x[0];
  // define thefunction, depending on the location of the xVal
  // otherwise return an error
  if( xVal >= mCompDomain[0][0] && xVal < mIntDomain[0][0] )
  {
    Double delta = mIntDomain[0][0]-mCompDomain[0][0];
    return mC*pow((mIntDomain[0][0]-xVal)/delta,3.0);
  }
  else if( xVal >= mIntDomain[0][0] && xVal <= mIntDomain[0][1] )
  {
    return 0.0;
  }
  else if( xVal > mIntDomain[0][1] && xVal < mCompDomain[0][1] )
  {
    Double delta = mCompDomain[0][1]-mIntDomain[0][1];
    return mC*pow((xVal-mIntDomain[0][1])/delta,3.0);
  }
  else
  {
    cerr << "ERROR: Evaluation point out of range!" << endl;
    exit(1);
  }
}

// the derivative of sigma1
Double PMLCubic2D::sigma1Dx( DoubleVector x, Double omega ) const
{
  // check if the input is good
  assert( x.size()==2 );

  // optain the x value
  Double xVal = x[0];
  // define thefunction, depending on the location of the xVal
  // otherwise return an error
  if( xVal >= mCompDomain[0][0] && xVal < mIntDomain[0][0] )
  {
    Double delta = mIntDomain[0][0]-mCompDomain[0][0];
    return -3.0*mC*pow(mIntDomain[0][0]-xVal,2.0)/pow(delta,3.0);
  }
  else if( xVal >= mIntDomain[0][0] && xVal <= mIntDomain[0][1] )
    return 0.0;
  else if( xVal > mIntDomain[0][1] && xVal < mCompDomain[0][1] )
  {
    Double delta = mCompDomain[0][1]-mIntDomain[0][1];
    return 3.0*mC*pow(xVal-mIntDomain[0][1],2.0)/pow(delta,3.0);
  }
  else
  {
    cerr << "ERROR: Evaluation point out of range!" << endl;
    exit(1);
  }
}

// the auxiliary function sigma2 in the PML 
Double PMLCubic2D::sigma2( DoubleVector x, Double omega ) const
{
  // check if the input is good
  assert( x.size()==2 );

  // optain the y value
  Double xVal = x[1];

  // define thefunction, depending on the location of the xVal
  // otherwise return an error
  if( xVal >= mCompDomain[1][0] && xVal < mIntDomain[1][0] )
  {
    Double delta = mIntDomain[1][0]-mCompDomain[1][0];
    return mC*pow((mIntDomain[1][0]-xVal)/delta,3.0);
  }
  else if( xVal >= mIntDomain[1][0] && xVal <= mIntDomain[1][1] )
    return 0.0;
  else if( xVal > mIntDomain[1][1] && xVal < mCompDomain[1][1] )
  {
    Double delta = mCompDomain[1][1]-mIntDomain[1][1];
    return mC*pow((xVal-mIntDomain[1][1])/delta,3.0);
  }
  else
  {
    cerr << "ERROR: Evaluation point out of range!" << endl;
    exit(1);
  }
}

// the derivative of sigma2
Double PMLCubic2D::sigma2Dy( DoubleVector x, Double omega ) const
{
  // check if the input is good
  assert( x.size()==2 );

  // optain the y value
  Double xVal = x[1];

  // define thefunction, depending on the location of the xVal
  // otherwise return an error
  if( xVal >= mCompDomain[1][0] && xVal < mIntDomain[1][0] )
  {
    Double delta = mIntDomain[1][0]-mCompDomain[1][0];
    return -3.0*mC*pow(mIntDomain[1][0]-xVal,2.0)/pow(delta,3.0);
  }
  else if( xVal >= mIntDomain[1][0] && xVal <= mIntDomain[1][1] )
    return 0.0;
  else if( xVal > mIntDomain[1][1] && xVal < mCompDomain[1][1] )
  {
    Double delta = mCompDomain[1][1]-mIntDomain[1][1];
    return 3.0*mC*pow(xVal-mIntDomain[1][1],2.0)/pow(delta,3.0);
  }
  else
  {
    cerr << "ERROR: Evaluation point out of range!" << endl;
    exit(1);
  }
}

