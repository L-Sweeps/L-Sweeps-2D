//include "Problem.h"
#include <iostream>

#include "DD.h"
#include "RHS.h"
#include "GlobalIterativeSolver.h"
#include "Output.h"
#include "Input.h"
#include "OutputVTK.h"
#include <unistd.h>
#include <sys/time.h>

using namespace std;

// The standard constructor, generates an empty object
Problem::Problem(): mDD(NULL), mRHS(NULL), mSolver(NULL), mOutput(NULL)
{

}

// The copy constructor
Problem::Problem( const Problem &other ): mDD(other.mDD), mRHS(other.mRHS), mSolver( other.mSolver), mOutput(other.mOutput)
{

}

// The destructor, frees all memory
Problem::~Problem()
{
  if( mDD!=NULL )
    delete mDD;
  if( mRHS!=NULL )
    delete mRHS;
  if( mSolver!=NULL )
    delete mSolver;
  if( mOutput!=NULL )
    delete mOutput;
}

// The assignment operator
Problem &Problem::operator=( const Problem &other )
{
  // only assign if other is different from this object
  if( this != &other )
  {
    // for each member variable, free the stored memory if there is any
    // and copy the new values from toher

    if( mDD !=NULL )
      delete mDD;
    mDD = other.mDD->copy();

    if( mRHS !=NULL )
      delete mRHS;
    mRHS = other.mRHS->copy();

    if( mSolver !=NULL )
      delete mSolver;
    mSolver = other.mSolver->copy();

    if( mOutput !=NULL )
      delete mOutput;
    mOutput = other.mOutput->copy();
  }

  return *this;
}

// The function that reads a problem from an Input class
void Problem::read( const Input &input, UInt probID )
{
  // for each member variable, free memory of possibly saved values, then read each member from input
  if( mDD!=NULL )
    delete mDD;
  input.readProblemDD( mDD, probID );

  if( mRHS!=NULL )
    delete mRHS;
  input.readProblemRHS( mDD, mRHS, probID );

  if( mSolver!=NULL )
    delete mSolver;
  input.readProblemSolver( mSolver, probID );
  mSolver->setParent( mDD );

  if( mOutput!=NULL )
    delete mOutput;
  input.readProblemOutput( mOutput, probID );
  if( mOutput != NULL )
    mOutput->setDD( mDD );
}

// The function that solves the problem
void Problem::solve() const
{
  // get the current rank, this is used for output to screen
  // only rank 0 outputs to screen
  UInt rank;
  MPI_Comm_rank(mDD->getComm(), (int*)&rank);

  // allocate variables for timing
  long start, end;
  struct timeval timecheck;
  double factorize, precond, solve;
  // if rank is 0, start time measurement and print to screen
  if( rank==0 )
  {
    cout << "Facorize problem..." << endl;
    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
  }
  // factorize the problem
  mDD->factorize();
  // if rank is 0, measure time and print to screen
  if( rank==0 )
  {
    gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    factorize = (end-start)/1000.0;
    
    cout << "Assemble RHS..." << endl;
  }
  // assemble the right hand sied
  vector<ComplexVector> rhs;
  mDD->assembleRHS( mRHS, rhs );
  // if rank is 0, solve the problem and measure time
  if( rank==0 )
  {
    cout << "Solving..." << endl;
    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;
  }
  vector<ComplexVector> res;
  mSolver->solve(rhs, res);
  // if rank is 0, print out the info
  if( rank==0 )
  {
    gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    solve = (end-start)/1000.0;
    cout << "Time to factorize system: " << factorize << endl;
    cout << "Time to run the GMRES: " << solve << endl;
  }

  // now write out the output file if Output is not a NULL vector
  if( rank==0 )
  {
    cout << "Writing output file..." << endl;
  }
  if( mOutput != NULL )
  {
    mOutput->setVals(res);
    mOutput->print();
  }

}
