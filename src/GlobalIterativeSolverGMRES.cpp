#include <iostream>
#include "GlobalIterativeSolverGMRES.h"
#include <mpi.h>

using namespace::std;

// the standard constructor
// we just set some standard values for tolerance, maximum number of iterations and verboseness
GlobalIterativeSolverGMRES::GlobalIterativeSolverGMRES(): GlobalIterativeSolver(), mTol(1.0e-6), mMaxIter(40), mVerbose(true), mLog(0)
{

}

// the copy constructor
GlobalIterativeSolverGMRES::GlobalIterativeSolverGMRES( const GlobalIterativeSolverGMRES &other ): GlobalIterativeSolver(other), mTol(other.mTol), mMaxIter(other.mMaxIter), mVerbose(other.mVerbose), mLog(other.mLog)
{

}

// the assignment operator
GlobalIterativeSolverGMRES &GlobalIterativeSolverGMRES::operator=( const GlobalIterativeSolverGMRES &other )
{
  if( &other != this )
  {
    mParent = other.mParent;
    mTol = other.mTol;
    mMaxIter = other.mMaxIter;
    mVerbose = other.mVerbose;
    mLog = other.mLog;
  }
  return *this;
}

// the function to set the tolrance
void GlobalIterativeSolverGMRES::setTol( Double tol )
{
  mTol = tol;
}

// the function to set the maximum number of iterations
void GlobalIterativeSolverGMRES::setMaxIter( UInt maxIter )
{
  mMaxIter = maxIter;
}

// the function to set the verboseness
void GlobalIterativeSolverGMRES::setVerbose( bool verbose )
{
  mVerbose = verbose;
}

// the function to return the log of the GMRES method
void GlobalIterativeSolverGMRES::getLog( vector<Double> &res ) const
{
  res = mLog;
}

// the copy function to return a GlobalIterativeSolver pointer to this object
GlobalIterativeSolver *GlobalIterativeSolverGMRES::copy() const
{
  return new GlobalIterativeSolverGMRES(*this);
}

// the function that realizes the solve
void GlobalIterativeSolverGMRES::solve(const vector<ComplexVector> &rhs, vector<ComplexVector> &res) const
{
  // get the communicator and the rank
  UInt rank;
  MPI_Comm comm = mParent->getComm();
  MPI_Comm_rank( comm, (int*)&rank );

  // get the right hand side for the preconditioned GMRES by
  // applying the preconditioner to the right hand side
  vector<ComplexVector> rhsPrec;
  mParent->applyPrecond( rhs, rhsPrec );

  // start with a zero initial guess
  res = rhsPrec;

  vector<ComplexVector> r, rTemp;
  mParent->applySystemMatrix(res, rTemp);
  mParent->applyPrecond( rTemp, r );
  for( UInt i=0; i<r.size(); i++ )
    r[i] = rhsPrec[i]-r[i];

  // compute the norm of the first residual (which is just the rhs)
  Double rhs_normlocsqrd = 0.0;
  Double r_normlocsqrd = 0.0;
  UInt localSize = 0;
  // compute the local size of the residual and the local norms
  for( UInt k=0; k<rhsPrec.size(); k++ )
  {
    rhs_normlocsqrd += pow(rhsPrec[k].norm(), 2.0);
    r_normlocsqrd += pow(r[k].norm(), 2.0);
    localSize += rhsPrec[k].size();
  }
  MPI_Barrier(comm);

  // communicate between all processes to compute the norm
  Double rhs_norm, r_norm;
  MPI_Allreduce( &rhs_normlocsqrd, &rhs_norm, 1, MPI_DOUBLE, MPI_SUM, comm ); 
  MPI_Allreduce( &r_normlocsqrd, &r_norm, 1, MPI_DOUBLE, MPI_SUM, comm ); 
  MPI_Barrier(comm);
  rhs_norm = pow(rhs_norm,0.5);
  r_norm = pow(r_norm,0.5);
  
  // save the residual information
  mLog = vector<Double>(1);
  mLog[0] = r_norm/rhs_norm;
  if( mVerbose && rank==0 )
    cout << "iteration 0: " << mLog[mLog.size()-1] << endl;

  // if this residual is already small enough, we are done
  if( mLog[0] < mTol || mMaxIter==0 )
    return;

  // allocate memory for the iterations
  // and initialize the data
  ComplexVector sn(0);
  ComplexVector cs(0);
  vector<ComplexMatrix> Q(r.size());
  for( UInt k=0; k<r.size(); k++ )
  {
    Q[k] = ComplexMatrix(r[k].size(),1);
    for( UInt l=0; l<r[k].size(); l++ )
      Q[k](l,0) = r[k][l]/r_norm; 
  }
  vector<Complex> beta(1);
  beta[0] = r_norm;
  ComplexMatrix H(0,0);
  UInt lastk = 0;

  // start the iterations (we do at most mMaxIter iterations)
  for( UInt k=0; k<mMaxIter; k++ )
  {
    // apply the arnoldi step
    ComplexVector h;
    vector<ComplexVector> q;
    arnoldi(comm,Q,k,h,q);

    // apply the Givens rotations
    apply_givens_rotation(h,cs,sn,k);
    MPI_Barrier(comm);

    // update H
    if( k==0 )
    {
      H = ComplexMatrix(h.size(),1);
      for( UInt i=0; i<h.size(); i++ )
        H(i,0) = h[i];
    }
    else
    {
      H.conservativeResize(H.rows()+1,H.cols()+1);
      for( UInt i=0; i<H.cols()-1; i++ )
        H(H.rows()-1,i) = 0.0;
      for( UInt i=0; i<H.rows(); i++ )
        H(i,H.cols()-1) = h[i];
    }

    // update Q
    for( UInt i=0; i<Q.size(); i++ )
    {
     Q[i].conservativeResize(Q[i].rows(), Q[i].cols()+1);
     for( UInt j=0; j<Q[i].rows(); j++ )
       Q[i](j,Q[i].cols()-1) = q[i][j];
    }
    MPI_Barrier(comm);

    // update beta
    beta.resize( beta.size()+1 );
    beta[beta.size()-1] = -sn[k]*beta[k];
    beta[k] *=cs[k];
    MPI_Barrier(comm);

    // save the residual info
    mLog.push_back( abs(beta[beta.size()-1]) / rhs_norm );
    MPI_Barrier(comm);

    // if verbose is true, print the residual and iteration count
    if( mVerbose && rank==0 )
      cout << "iteration " << k+1 << ": " << mLog[mLog.size()-1] << endl;

    // if the resdiual is small enough, stop
    if( mLog[mLog.size()-1] <= mTol || k==mMaxIter-1 )
    {
      lastk = k;
      break;
    }
  }
  MPI_Barrier(comm);

  // reconstruct the solution
  UInt k = lastk;
  ComplexVector rhsTemp(k+1);
  for( UInt i=0; i<k+1; i++ )
    rhsTemp[i] = beta[i];

  ComplexVector y = H.block(0,0,k+1,k+1).colPivHouseholderQr().solve(rhsTemp);
  for( UInt i=0; i<res.size(); i++ )
    res[i] += Q[i].block(0,0,Q[i].rows(),k+1)*y;
  MPI_Barrier(comm);
}

// the function that converts a vector of vectors over subdomains to a global vector
// no memory is allocated, it is assumed that the masks of the given vectors is correct
void GlobalIterativeSolverGMRES::vectorize( MPI_Comm comm, const vector<ComplexVector> &x, ComplexVector &res ) const
{
  // the current index of the global vector
  UInt curr_ind=0;
  // loop over all small vectors and append
  for( UInt k=0; k<x.size(); k++ )
    for(UInt l=0; l<x[k].size(); l++ )
    {
      res[curr_ind] = x[k][l];
      curr_ind++;
    }
  MPI_Barrier(comm);
}

// the function that converts a global vector to a vector of vectors
// no memory is allocated, it is assumed that the masks of the given vectors is correct
void GlobalIterativeSolverGMRES::devectorize( MPI_Comm comm, const ComplexVector &x, vector<ComplexVector> &res ) const
{
  // the current index of the global vector
  UInt curr_ind=0;
  // loop over all small vectors and set them up
  for( UInt k=0; k<res.size(); k++ )
    for(UInt l=0; l<res[k].size(); l++ )
    {
      res[k][l] = x[curr_ind];
      curr_ind++;
    }
  MPI_Barrier(comm);
}

// the function that realizes the arnoldi step
void GlobalIterativeSolverGMRES::arnoldi( MPI_Comm comm, const vector<ComplexMatrix> &Q, UInt k, ComplexVector &h, vector<ComplexVector> &q ) const
{
  UInt rank;
  MPI_Comm_rank( comm, (int*)&rank );
  MPI_Barrier(comm);
  // get the k-th colum of Q
  vector<ComplexVector> QTemp(Q.size());
  for( UInt i=0; i<Q.size(); i++ )
  {
    QTemp[i] = ComplexVector(Q[i].rows());
    for( UInt j=0; j<Q[i].rows(); j++ )
      QTemp[i][j] = Q[i](j,k);
  }
  // apply the system matrix (and the preconditioner) to the k-th column 
  mParent->applySystemMatrix(QTemp, q);
  QTemp = q;
  MPI_Barrier(comm);
  mParent->applyPrecond(QTemp, q );

  // compute the size of the global vector to allocate memory
  UInt size = 0;
  for( UInt i=0; i<q.size(); i++ )
    size += q[i].size();
  // allocate memory and convert the vector of vectors to a global vector
  ComplexVector qVec(size);
  vectorize(comm, q,qVec);
  MPI_Barrier(comm);

  // allocate memory and start the loop
  h = ComplexVector(k+2);
  for( UInt i=0; i<k+1; i++ )
  {
    MPI_Barrier(comm);
    // get the i-th column of Q
    for( UInt j=0; j<Q.size(); j++ )
    {
      for( UInt l=0; l<Q[j].rows(); l++ )
        QTemp[j][l] = Q[j](l,i);
    }
    // convert the i-th column to a global vector
    ComplexVector QTempVec(size);
    vectorize(comm,QTemp,QTempVec);

    // update h and q
    Complex hiiloc = qVec.dot(QTempVec);
    Double hiiloc_re = real(hiiloc);
    Double hiiloc_im = imag(hiiloc);
    Double hii_re, hii_im;
    MPI_Barrier(comm);
    MPI_Allreduce(&hiiloc_re, &hii_re, 1, MPI_DOUBLE, MPI_SUM, comm ); 
    MPI_Allreduce(&hiiloc_im, &hii_im, 1, MPI_DOUBLE, MPI_SUM, comm );
    MPI_Barrier(comm);
    Complex im(0.0,1.0);
    h[i] = hii_re+im*hii_im;
    qVec -= h[i]*QTempVec;
  }
  MPI_Barrier(comm);
  // finish updateing h and q
  Double q_norm_loc_sqrd = pow(qVec.norm(),2.0);
  MPI_Barrier(comm);
  Double q_norm;
  MPI_Allreduce(&q_norm_loc_sqrd,&q_norm,1,MPI_DOUBLE,MPI_SUM,comm);
  MPI_Barrier(comm);
  q_norm = pow(q_norm,0.5);
  MPI_Barrier(comm);
  h[k+1] = q_norm;
  qVec/=h[k+1];
    

  // convert the global vector of q back to the vector of vectors defined over subdomains
  devectorize(comm,qVec,q);
  MPI_Barrier(comm);
}

// the function to apply the givens rotiation
void GlobalIterativeSolverGMRES::apply_givens_rotation(ComplexVector &h, ComplexVector &cs, ComplexVector &sn, UInt k ) const
{
  for( UInt i=1; i<k; i++ )
  {
    Complex temp = cs[i-1]*h[i-1] + sn[i-1]*h[i];
    h[i] = -sn[i-1]*h[i-1] + cs[i-1]*h[i];
    h[i-1] = temp;
  } 

  ComplexVector temp = givens_rotation( h[k], h[k+1] );
  h[k] = temp[0]*h[k]+temp[1]*h[k+1];
  h[k+1] = 0.0;

  cs.conservativeResize(cs.size()+1);
  cs[cs.size()-1] = temp[0];
  sn.conservativeResize(sn.size()+1);
  sn[sn.size()-1] = temp[1];
}

// the functions that realizes one givens rotation
ComplexVector GlobalIterativeSolverGMRES::givens_rotation( Complex v1, Complex v2 ) const
{
  ComplexVector res(2);
  if(v1==(0.0,0.0))
  {
    res[0] = 0.0;
    res[1] = 1.0;
  }
  else
  {
    Complex t = pow( pow(v1,2.0)+pow(v2,2.0), 0.5 );
    res[0] = abs(v1) / t;
    res[1] = res[0]*v2 / v1;
  }
  return res;
}
