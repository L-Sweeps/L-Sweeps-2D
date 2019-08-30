#include<iostream>
#include <assert.h>

#include "FiniteDifferences.h"
#include "typedef.h"
#include "utils.h"

using namespace::std;

// implements the Fornberg Algorithm
DoubleVector FDweights( Int z, IntVector x, Int m)
{
  Int n = x.rows()-1;
  Double c1 = 1.0;
  Double c4 = x[0]-z;

  DoubleVector c( (n+1)*(m+1) );
  for( UInt k=0; k<c.size(); k++ )
    c[k] = 0.0;
  c[0] = 1.0;
  for( Int i=0; i<n; i++ )
  {
    Int mn = min(i+1,m);
    Double c2 = 1.0;
    Double c5 = c4;
    c4 = x[i+1]-z;

    for( Int j=0; j<=i; j++ )
    {
      Double c3 = x[i+1]-x[j];
      c2 *= c3;
      if( j==i )
      {
        for( Int k=mn-1; k>=0; k-- )
        {
          c[(k+1)*(n+1)+i+1] = c1*((k+1)*c[k*(n+1)+i]-c5*c[(k+1)*(n+1)+i])/c2;
        }
        c[i+1] = -c1*c5*c[i]/c2;
      }
      for( Int k=mn-1; k>=0; k-- )
      {
        c[(k+1)*(n+1)+j] = (c4*c[(k+1)*(n+1)+j]-(k+1)*c[k*(n+1)+j])/c3;
      }
      c[j] = c4*c[j]/c3;
    }
    c1=c2;
  }

  DoubleVector cwei( n+1 );
  for( Int i=0; i<n+1; i++ )
    cwei[i] = c[m*(n+1)+i];

  return cwei;
}

// assemble the finite difference matrix to approximate derivatives in 1D
DoubleSparseMatrix DerivativeDifferenceMatrix1D( UInt nx, Double h, UInt approxOrder, UInt derOrder )
{
  // there have to be at least order+1 points
  assert( nx>=approxOrder+1 );

  // allocate memory to save the entry information of the sparse matrix
  DoubleVector vals( (approxOrder+1)*nx );
  UIntVector rowInds( (approxOrder+1)*nx );
  UIntVector colInds( (approxOrder+1)*nx );

  // find the center of the stencil
  Int center = approxOrder/2;

  // collect the evaluation points of the stencil
  IntVector intVec(approxOrder+1);
  for( UInt i=0; i<approxOrder+1; i++ )
    intVec[i] = i;

  // assemble the stencils on the left end of the boundary
  for( Int i=0; i<center; i++ )
  {
    // compute the local stencil
    DoubleVector locCoeffs = FDweights( i, intVec, derOrder );

    // add it to the entries
    for( Int k=0; k<approxOrder+1; k++ )
    {
      vals[i*(approxOrder+1)+k] = locCoeffs[k];
      rowInds[i*(approxOrder+1)+k] = i;
      colInds[i*(approxOrder+1)+k] = k;
    }
  }

  // assemble the middle stencils
  // precompute the stencil (it is the same for all these points)
  DoubleVector locCoeffs = FDweights( center, intVec, derOrder );
  // add the entries
  for( Int i=max(0,center); i<nx-max(center,Int(approxOrder-center)); i++ )
  {
    for( Int k=0; k<approxOrder+1; k++ )
    {
      vals[i*(approxOrder+1)+k] = locCoeffs[k];
      rowInds[i*(approxOrder+1)+k] = i;
      colInds[i*(approxOrder+1)+k] = i-max(center,0)+k;
    }
  }

  // assemble the stencils on the right end of the boundary
  for( Int i=nx-max(center,Int(approxOrder-center)); i<nx; i++ )
  {
    // compute the local stencil
    DoubleVector locCoeffs = FDweights( nx-1-i+approxOrder, intVec, derOrder );

    // add it to the entries
    for( Int k=0; k<approxOrder+1; k++ )
    {
      vals[i*(approxOrder+1)+k] = locCoeffs[k];
      rowInds[i*(approxOrder+1)+k] = i;
      colInds[i*(approxOrder+1)+k] = nx-(approxOrder+1)+k;
    }
  }
  // add the scaling
  vals /= pow(h,derOrder);

  return generateSparseMatrix( nx, nx, rowInds, colInds, vals );
}
