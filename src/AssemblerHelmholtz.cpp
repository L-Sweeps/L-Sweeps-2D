#include "AssemblerHelmholtz.h"
#include "Function.h"
#include <iostream>

using namespace::std;

// sets up an empty object, makes sure that mPML is set to NULL
AssemblerHelmholtz::AssemblerHelmholtz(): mPML(NULL), mM(NULL)
{

}

// frees the memory in mPML before destruction
AssemblerHelmholtz::~AssemblerHelmholtz()
{
  if( mPML != NULL )
    delete mPML;
  if( mM != NULL )
    delete mM;
}

// copies all member variables, especially mPML which is properly copied
AssemblerHelmholtz::AssemblerHelmholtz( const AssemblerHelmholtz &other ): Assembler( other ), mOmega( other.mOmega )
{
  // initialize PML and squared slowness
  mPML = NULL;
  mM = NULL;

  // if there is a PML, copy it
  if( other.mPML != NULL )
  {
    mPML = other.mPML->copy();
    mPML->setParent(this);
  }
  // if there is a squared slowness, copy it
  if( other.mM != NULL )
  {
    mM = other.mM->copy();
  }
}

// copies all information from other to this object
AssemblerHelmholtz &AssemblerHelmholtz::operator=( const AssemblerHelmholtz &other )
{
  // only copy if the two objects are not the same
  if( this != &other )
  {
    // assign as an Assembler object
    (Assembler &)(*this) = (const Assembler &)other;

    // free memory if pml and squared slowness are set in this object
    if( mPML != NULL )
      delete[] mPML;
    if( mM != NULL )
      delete[] mM;

    // initialize pml and squared slowness
    mPML = NULL;
    mM = NULL;

    // copy the pml from ohter
    mPML = other.mPML->copy();
    // properly set this object as the parent object of the PML
    mPML->setParent(this);

    // if there is a PML, copy it
    if( other.mPML != NULL )
    {
      mPML = other.mPML->copy();
      mPML->setParent(this);
    }
    // if there is a squared slowness, copy it
    if( other.mM != NULL )
    {
      mM = other.mM->copy();
    }
    
    // set the frequency omega
    mOmega = other.mOmega;
  }
  return *this;
}

// sets the PML
void AssemblerHelmholtz::setPML( PML *pml )
{
  // if mPML was previously set, free the memory
  if( mPML != NULL )
    delete mPML;
  // set the PML
  mPML = pml->copy();

  // oproperly set this object as the parent object of the PML
  mPML->setParent(this);
}

// sets the function m
void AssemblerHelmholtz::setM( Function* m )
{
  // if the squared slowness is set, delete the old one
  if( mM != NULL )
    delete mM;
  // set the squared slowness m
  mM = m->copy();
}

// sets the frequency omega
void AssemblerHelmholtz::setOmega( Double omega )
{
  mOmega = omega;
}
