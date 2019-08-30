/// \file utils.h
/// A file containing some utility functions.

#ifndef UTILS_H
#define UTILS_H

#include "typedef.h"

using namespace std::placeholders;
using namespace::std;

Double const pi=M_PI;

//! Generates a \type DoubleSparseMatrix from index information and non-zero entries.
/*! @param rows The number of rows of the matrix
 *  @param cols The number of columns of the matrix
 *  @param rowInds The vector containing the row indices, same length as number of non-zeros
 *  @param colInds The vector containing the column indices, same length as number of non-zeros
 *  @param vals The vector containing the matrix entries
 *  @return The generated \type DoubleSparseMatrix 
 */
DoubleSparseMatrix generateSparseMatrix( UInt rows, UInt cols, UIntVector rowInds, UIntVector colInds, DoubleVector vals );

//! Generates a \type ComplexSparseMatrix from index information and non-zero entries.
/*! @param rows The number of rows of the matrix
 *  @param cols The number of columns of the matrix
 *  @param rowInds The vector containing the row indices, same length as number of non-zeros
 *  @param colInds The vector containing the column indices, same length as number of non-zeros
 *  @param vals The vector containing the matrix entries
 *  @return The generated \type ComplexSparseMatrix 
 */
ComplexSparseMatrix generateSparseMatrix( UInt rows, UInt cols, UIntVector rowInds, UIntVector colInds, ComplexVector vals );


//! Reads in a \DoubleSparseMatrix from a text file.
/*! @param filename The filename to be read from
 *  @return The generated \type DoubleSparseMatrix 
 *
 * The file has to contain a block of the following structure:
 * ~~~~~~~~~~~~~~~~  
 * SparseMatrix
 * number of rows
 * number of cols
 * number of non-zero entries
 * rowIndex1
 * rowIndex2
 * ...
 * rowIndexN
 * colIndex1
 * colIndex2
 * ...
 * colIndedN
 * val1
 * val2
 * ...
 * valN
 * ~~~~~~~~~~~~~~~~  
 */
DoubleSparseMatrix readSparseMatrix( string filename );

//! Computes the kronecker product of two sparse matrices
/*!
 * Computes \f$ A \otimes B\f$.
 * \param A The matrix \f$A\f$
 * \param B The matrix \f$B\f$
 * \return The result of the product
 */
ComplexSparseMatrix kron( const ComplexSparseMatrix &A, const ComplexSparseMatrix &B );

//! Computes the kronecker product of two sparse matrices
/*!
 * Computes \f$ a \otimes b\f$.
 * \param a The matrix \f$a\f$
 * \param b The matrix \f$b\f$
 * \return The result of the product
 */
ComplexVector kron( const ComplexVector &a, const ComplexVector &b );

//! Slices a matrix for given row and column indices
/*!
 * \param A The matrix to be sliced
 * \param rowInds The row indices to use in the slice
 * \param colInds The column indices to use in the slice
 * \retrn The sliced matrix
 */
ComplexSparseMatrix slice( const ComplexSparseMatrix &A, UIntVector rowInds, UIntVector colInds );

//! Sorts a vector and returns the transformation of the sorting (in indices of the orignal vector)
/*!
 * \param x The vector to be sorted
 * \para res The transformation, given as a vector of original indices of \p x
 */
void getSortedIndices( const DoubleVector &x, UIntVector &res );

//! Generates a Julia file that can be run to produce a plot using physical DOF points and a vector of values
/*!
 * \param pts The physical DOF points
 * \param vals The values to be ploted, it uses the real value
 * \param filename The .jl file where the plotting algorithm should be stored
 * \param filename_output The filename of the .png that should be produced by the .jl file.
 */
void generateJulia2DRealPlot( vector<DoubleVector> pts, ComplexVector vals, string filename, string output_filename );

#endif
