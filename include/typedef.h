/// \file typedef.h
///   A file containing all basic typedefs.
///

#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <complex>
#include <assert.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace::std;


// Integer typedefs


//! Defines the integer type used throughout the code. Just overwrites the redefines C++ int.
typedef int Int;

//! Defines a vector type of \type Int using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::Matrix<Int,Eigen::Dynamic,1> IntVector; 

//! Defines a Matrix type of \type Int using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::Matrix<Int,Eigen::Dynamic,Eigen::Dynamic> IntMatrix; 



// Unsigned Integer typedefs


//! Defines the unsigned integer type used throughout the code. Just redefines the standard C++ unsigned int.
typedef unsigned int UInt;

//! Defines a vector type of \type UInt using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::Matrix<UInt,Eigen::Dynamic,1> UIntVector; 

//! Defines a matrix type of \type UInt using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::Matrix<UInt,Eigen::Dynamic,Eigen::Dynamic> UIntMatrix; 



// Floating Point typedefs


//! Defines the type for floating pont numbers used throughout the code. Just redefines the standard C++ double.
typedef double Double;

//! Defines a vector type of \type Double using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::Matrix<Double,Eigen::Dynamic,1> DoubleVector; 

//! Defines a matrix type of \type Double using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::Matrix<Double,Eigen::Dynamic,Eigen::Dynamic> DoubleMatrix; 

//! Defines a sparse matrix type of \type Double using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::SparseMatrix<Double> DoubleSparseMatrix; 


// Complex number typedefs


/// Defines the type for complex numbers used throughout the code. Redefines the standard C++ complex template, for details see [here](link https://en.cppreference.com/w/cpp/numeric/complex).
typedef std::complex<Double> Complex; 

//! Defines a vector type of \type Complex using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::Matrix<Complex,Eigen::Dynamic,1> ComplexVector; 

//! Defines a matrix type of \type Complex using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic> ComplexMatrix; 

//! Defines a sparse matrix type of \type Double using [Eigen](https://eigen.tuxfamily.org/dox/).
typedef Eigen::SparseMatrix<Complex> ComplexSparseMatrix;

//! Defines a type for a function pointer
typedef Double (*FunctionPtr)(DoubleVector);

#endif
