#include "GlobalIterativeSolver.h"

// the standard constructor
GlobalIterativeSolver::GlobalIterativeSolver(): mParent(NULL)
{

}

// the copy constructor
GlobalIterativeSolver::GlobalIterativeSolver( const GlobalIterativeSolver &other ): mParent(other.mParent)
{

}

// the assignment operator
GlobalIterativeSolver &GlobalIterativeSolver::operator=( const GlobalIterativeSolver &other )
{
  if( &other != this )
  {
    mParent = other.mParent;
  }
  return *this;
}

// the function to set the parent DD
void GlobalIterativeSolver::setParent( const DD* parent )
{
  mParent = parent;
}
