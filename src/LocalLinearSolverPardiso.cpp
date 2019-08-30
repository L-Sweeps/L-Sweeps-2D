#include "LocalLinearSolverPardiso.h"
#include <iostream>
#include <float.h>

using namespace::Eigen;
using namespace::std;

//extern "C" void pardisoinit (void   *, int *,   int *, int *, Double *, int *);
//extern "C" void pardiso     (void   *, int *,   int *, int *,    int *, int *, 
//                  Complex *, int *,   int *, int *,   int *, int *,
//                  int *, Complex *, Complex *, int *, Double *);


extern "C" void pardisoinit (void   *, int *,   int *, int *, Double *, int *);
extern "C" void pardiso     (void   *, int *,   int *, int *,    int *, int *, 
                  Complex *, int *,   int *, int *,   int *, int *,
                  int *, Complex *, Complex *, int *, Double *);
extern "C" void pardiso_chkmatrix_z  (int *, int *, Complex *, int *, int *, int *);
extern "C" void pardiso_chkvec_z     (int *, int *, Complex *, int *);
extern "C" void pardiso_printstats_z (int *, int *, Complex *, int *, int *, int *,
                           Complex *, int *);

LocalLinearSolverPardiso::LocalLinearSolverPardiso(): mA(NULL), mIsFactorized(false), mIA(NULL), mJA(NULL), mHasMatrix(false)
{
  int mtype = 13;
  int solver = 0;
  int error = 0;

  int num_procs;
  char *var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
    sscanf( var, "%d", &num_procs );
  else {
        cerr << "Set environment OMP_NUM_THREADS" << endl;
        exit(1);
  }
  for(UInt i=0;i<64;i++)
    mIParm[i] = 0;
  mIParm[0] = 1;
  mIParm[1] = 2;
  mIParm[2]  = num_procs;
  mIParm[9]  = 13;
  mIParm[10]  = 1; 

  pardisoinit(mPt, &mtype, &solver, mIParm, mDParm, &error);

  if( error != 0 )
  {
    if( error == -10 )
      cerr << "No license file found for Pardiso" << endl;
    else if( error == -11 )
      cerr << "Pardiso license is expired" << endl;
    else if( error == -12 )
      cerr << "Wrong username or hostname in Pardiso licens" << endl;
    
    exit(1);
  }
  clear();
}

// implementation of the (private) copy constructor
LocalLinearSolverPardiso::LocalLinearSolverPardiso( const LocalLinearSolverPardiso &other ): mIsFactorized(false), mIA(NULL), mJA(NULL), mHasMatrix(false), mA(NULL)
{
  int mtype = 13;
  int solver = 0;
  int error = 0;

  int num_procs;
  char *var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
    sscanf( var, "%d", &num_procs );
  else {
        cerr << "Set environment OMP_NUM_THREADS" << endl;
        exit(1);
  }
  for(UInt i=0;i<64;i++)
    mIParm[i] = 0;
  mIParm[0] = 1;
  mIParm[1] = 2;
  mIParm[2]  = num_procs;
  mIParm[9]  = 13;
  mIParm[10]  = 1; 

  pardisoinit(mPt, &mtype, &solver, mIParm, mDParm, &error);

  mA = other.mA;
  if( other.mHasMatrix )
    setMatrix( other.mA );
 
  if( other.mIsFactorized )
    factorize();
}
    
// implementation of the (private) copy constructor
LocalLinearSolverPardiso &LocalLinearSolverPardiso::operator=(const LocalLinearSolverPardiso &other)
{
  if( this != &other )
  {
    clear();

    int mtype = 13;
    int solver = 0;
    int error = 0;

    int num_procs;
    char *var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
      sscanf( var, "%d", &num_procs );
    else {
      cerr << "Set environment OMP_NUM_THREADS" << endl;
      exit(1);
    }
  for(UInt i=0;i<64;i++)
    mIParm[i] = 0;
  mIParm[0] = 1;
  mIParm[1] = 2;
  mIParm[2]  = num_procs;
  mIParm[9]  = 13;
  mIParm[10]  = 1; 

    pardisoinit(mPt, &mtype, &solver, mIParm, mDParm, &error);

    mHasMatrix = false;
    mIsFactorized = false;

    mA = other.mA;
    if( other.mHasMatrix )
      setMatrix( other.mA );

    if( other.mIsFactorized )
      factorize();
  }
  return *this;
}

// generate a copy of this object and return a point to it
LocalLinearSolver *LocalLinearSolverPardiso::copy() const
{
  return new LocalLinearSolverPardiso(*this);
}

// sets the matrix of the linear system, the ComplexSparse matrix
// has to be compressed before this function is called, this is inherited from Eigen
void LocalLinearSolverPardiso::setMatrix( ComplexSparseMatrix *A )
{
  // set the matrix
  mA = A;

  if( mIA != NULL )
    delete mIA;
  mIA = new int[mA->rows()+1];
  for( UInt i=0; i<mA->rows()+1; i++ )
    mIA[i] = mA->outerIndexPtr()[i]+1;

  if( mJA != NULL )
    delete mJA;
  mJA = new int[mIA[mA->rows()]];
  for( UInt i=0; i<mIA[mA->rows()]; i++ )
    mJA[i] = mA->innerIndexPtr()[i]+1;

  mIsFactorized=false;
  mHasMatrix=true;
}

// factorizes the system
void LocalLinearSolverPardiso::factorize()
{
  if( mIsFactorized )
  {
    return;
  }
  if( !mHasMatrix )
  {
    cerr << "Matrix has to be set before factorization!" << endl;
    exit(1);
  }

  int n = mA->rows();
  int phase = 11;
  int maxfct = 1;
  int mnum = 1;
  int msglvl = 0;
  int error = 0;
  int mtype = 13;
  int nrhs = 1;
  Complex ddum;
  int idum;

  pardiso( mPt, &maxfct, &mnum, &mtype, &phase, &n, mA->valuePtr(), mIA, mJA, &idum, &nrhs, mIParm, &msglvl, &ddum, &ddum, &error, mDParm );
  if( error!= 0 )
  {
    cerr << "ERROR: Pardiso factorization failed in symbolic factorization." << endl;
    exit(1);
  }

  phase = 22;
  pardiso( mPt, &maxfct, &mnum, &mtype, &phase, &n, mA->valuePtr(), mIA, mJA, &idum, &nrhs, mIParm, &msglvl, &ddum, &ddum, &error, mDParm);

  if( error!= 0 )
  {
    cerr << "ERROR: Pardiso factorization failed in numerical factorization." << endl;
    exit(1);
  }

  // set the flag that the system is factorized
  mIsFactorized = true;
}

// solves the linear system for the rhs b and returns the result in res
void LocalLinearSolverPardiso::solve(ComplexVector &b, ComplexVector &res) const
{
  // if the system is not factorized, do it
  if( !mIsFactorized )
  {
    cerr << "ERROR: The linear system has to be factorized before it can be solved!" << endl;
    exit(1);
  }

  int phase=33;
  mIParm[11]=1;

  int n = mA->rows();
  int maxfct = 1;
  int mnum = 1;
  int msglvl = 0;
  int error = 0;
  int mtype = 13;
  int nrhs = 1;
  Complex ddum;
  int idum;

  res = b;
  Double max = abs(b[0]);
  for( UInt i=0; i<b.size(); i++ )
  {
    Double val = abs( b[i] );
    if( val>max )
      max = val;
  }
  max += 1e-15;;
  b/=max;
  pardiso(mPt, &maxfct, &mnum, &mtype, &phase, &n, mA->valuePtr(), mIA, mJA, &idum, &nrhs, mIParm, &msglvl, b.data(), res.data(), &error, mDParm);
  
  res *= max;
  b *= max;

  if( error != 0 )
  {
    cerr << "ERROR: Pardiso backsubstitution failed! Error code: " << error << endl;
    exit(1);
  }
}

LocalLinearSolverPardiso::~LocalLinearSolverPardiso()
{
  clear();
}

void LocalLinearSolverPardiso::clear()
{
  if( mHasMatrix )
  {
    int n = mA->rows();

    int phase=-1;
    int maxfct = 1;
    int mnum = 1;
    int msglvl = 0;
    int error = 0;
    int mtype = 13;
    int nrhs = 1;
    Complex ddum;
    int idum;

    pardiso(mPt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, mIA, mJA, &idum, &nrhs, mIParm, &msglvl, &ddum, &ddum, &error, mDParm);

  }
  if( mIA != NULL )
    delete mIA;
  if( mJA != NULL )
    delete mJA;

  mIA=NULL;
  mJA=NULL;
}

