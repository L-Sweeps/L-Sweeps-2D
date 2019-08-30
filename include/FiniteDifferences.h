#ifndef FINITE_DIFFERENCES_H
#define FINITE_DIFFERENCES_H

#include "typedef.h"

//! Generates the coefficients for a finite difference approximation of a derivative
/*! @param z The the evaluation point.
  @param x The vetor of all evaluation poitns.
  @param m The order of the derivative.
  @return The \type DoubleVector of the coefficients.

  Uses the [Fornberg algorithm](https://pdfs.semanticscholar.org/8bf5/912bde884f6bd4cfb4991ba3d077cace94c0.pdf).
  */
DoubleVector FDweights( Int z, IntVector x, Int m);

//! Generates the finite difference matrix for the first derivative in 1D
/*!
 *  @param nx The number of evaluation point (size of the matrix), this includes the endpoints 
 *  @param h The step-size use for the discretiztaion
 *  @param approxOrder The order of the finite difference approximation
 *  @param derOrder The order of the derivative to be approximated 
 *  @return The finite difference matrix as a \type DoubleSparseMatrix
 *
 *  This function assumes a uniform discretization
 */
DoubleSparseMatrix DerivativeDifferenceMatrix1D( UInt nx, Double h, UInt approxOrder, UInt derOrder );

#endif
