#include "utils.h"
#include<fstream>
#include "typedef.h"
#include <iostream>

// generates a sparse matrix from given input
DoubleSparseMatrix generateSparseMatrix( UInt rows, UInt cols, UIntVector rowInds, UIntVector colInds, DoubleVector vals )
{
  // check if the input is consistent
  assert( rowInds.rows()==colInds.rows() );
  assert( colInds.rows()==vals.rows() );

  // reserve space for input for Eigen
  typedef Eigen::Triplet<Double> T;
  std::vector<T> tripletList;
  tripletList.reserve(rowInds.rows());

  // set up the Eigen input from the input vectors
  for( UInt k=0; k<rowInds.rows(); k++ )
  {
    tripletList.push_back( T(rowInds[k],colInds[k],vals[k]) );
  }

  // generate the Eigen matrix and set the entries
  DoubleSparseMatrix res( rows, cols );
  res.setFromTriplets( tripletList.begin(), tripletList.end() );

  return res;
}

// generates a sparse matrix from given input
ComplexSparseMatrix generateSparseMatrix( UInt rows, UInt cols, UIntVector rowInds, UIntVector colInds, ComplexVector vals )
{
  // check if the input is consistent
  assert( rowInds.rows()==colInds.rows() );
  assert( colInds.rows()==vals.rows() );

  // reserve space for input for Eigen
  typedef Eigen::Triplet<Complex> T;
  std::vector<T> tripletList;
  tripletList.reserve(rowInds.rows());

  // set up the Eigen input from the input vectors
  for( UInt k=0; k<rowInds.rows(); k++ )
  {
    tripletList.push_back( T(rowInds[k],colInds[k],vals[k]) );
  }

  // generate the Eigen matrix and set the entries
  ComplexSparseMatrix res( rows, cols );
  res.setFromTriplets( tripletList.begin(), tripletList.end() );

  return res;
}

// reads in a sparse matrix from a file
DoubleSparseMatrix readSparseMatrix( string filename )
{
  // open the file for reading
  ifstream file;
  file.open(filename);

  // define buffer variables for read in lines and a flag if a matrix block was found
  string line;
  bool foundMatrix = false;

  // define variables for the entry information of the sparse matrix
  UIntVector rowInds(0);
  UIntVector colInds(0);
  DoubleVector vals(0);
  UInt rows = 0;
  UInt cols = 0;

  // if the file is properly opened start reading in
  if( file.is_open() )
  {
    // until a matrix block is found or you reach the end of the file, read in the file line by line
    while( !foundMatrix && getline(file, line) )
    {
      // if the start of the matrix block is found, start reading in
      if( line=="SparseMatrix" )
      {
        // read in the number of rows
        getline(file,line);
        rows = stoi(line);
        // read in the number of cols
        getline(file,line);
        cols = stoi(line);

        // read in the number of non-zero entries
        getline(file,line);
        UInt N = stoi(line);

        // allocate the required memory
        rowInds = UIntVector(N);
        colInds = UIntVector(N);
        vals = DoubleVector(N);

        // read in all row indices 
        for( UInt i=0; i<N; i++ )
        {
          getline(file,line);
          UInt ind = stoi(line);
          rowInds[i] = ind; 
        }
        // read in all col indices 
        for( UInt i=0; i<N; i++ )
        {
          getline(file,line);
          UInt ind = stoi(line);
          colInds[i] = ind; 
        }
        // read in all values indices 
        for( UInt i=0; i<N; i++ )
        {
          getline(file,line);
          Double val = stod(line);
          vals[i] = val; 
        }
      }
    }
    //close the file and generate the matrix
    file.close();
    return generateSparseMatrix( rows, cols, rowInds, colInds, vals ); 
  }
  // if the file cannot be opened properly, stop and print an error
  else
  {
    cerr << "ERROR: File to be read from is not opened!" << endl;
    file.close();
    exit(1);
  }
}

// computes the kronecker product of two matrices
ComplexSparseMatrix kron( const ComplexSparseMatrix &A, const ComplexSparseMatrix &B )
{
  // compute the size of the resulting matrix
  UInt rows = A.rows()*B.rows();
  UInt cols = A.cols()*B.cols();

  // allocate memory for the row and column indices and the non-zero values
  UIntVector rowInds(A.nonZeros()*B.nonZeros());
  UIntVector colInds(A.nonZeros()*B.nonZeros());
  ComplexVector vals(A.nonZeros()*B.nonZeros());

  // loop over the non-zero entries and assigne the row and column indices and the values
  UInt index=0;
  for (UInt k=0; k<A.outerSize(); ++k)
  {
    for (ComplexSparseMatrix::InnerIterator itA(A,k); itA; ++itA)
    {
      for (UInt l=0; l<B.outerSize(); ++l)
      {
        for (ComplexSparseMatrix::InnerIterator itB(B,l); itB; ++itB)
        {
          rowInds[index] = itA.row()*B.rows()+itB.row();
          colInds[index] = itA.col()*B.cols()+itB.col();
          vals[index] = itA.value()*itB.value();
          index++;
        }
      }
    }
  }

  // generate the result matrix and return it
  return generateSparseMatrix( rows, cols, rowInds, colInds, vals );
}

// compute the kronecker product between two vectors
ComplexVector kron( const ComplexVector &a, const ComplexVector &b )
{
  // allocate the memory of the resulting vector
  ComplexVector res(a.rows()*b.rows());

  // loop over the two vectors and compute the result
  for (UInt k=0; k<a.rows(); ++k)
  {
    for (UInt l=0; l<b.rows(); ++l)
    {
      res[k*b.rows()+l] = a[k]*b[l];
    }
  }

  // return the result
  return res; 
}

// slices a matrix
ComplexSparseMatrix slice( const ComplexSparseMatrix &A, UIntVector rowInds, UIntVector colInds )
{
  // sort the matrix so that the slice appears in a rectangle in the upper left corner of the matrix
  //
  // compute the permutation matrices
  ComplexSparseMatrix P1(A.rows(),A.rows());
  ComplexSparseMatrix P2(A.cols(),A.cols());
  vector<bool> done1(A.rows());
  vector<bool> done2(A.cols());
  for( UInt i=0; i<A.rows(); i++ )
  {
    done1[i] = false;
  }
  for( UInt i=0; i<A.cols(); i++ )
  {
    done2[i] = false;
  }
  for( UInt i=0; i<rowInds.size(); i++ )
  {
    P1.insert(i,rowInds[i]) = 1.0;
    done1[rowInds[i]] = true;
  }
  for( UInt i=0; i<colInds.size(); i++ )
  {
    P2.insert(colInds[i],i) = 1.0;
    done2[colInds[i]] = true;
  }
  UInt index = rowInds.size();
  for( UInt i=0; i<A.rows(); i++ )
  {
    if(!done1[i])
    {
      P1.insert(index,i) = 1.0;
      index++;
    }
  }
  index = colInds.size();
  for( UInt i=0; i<A.cols(); i++ )
  {
    if(!done2[i])
    {
      P2.insert(i,index) = 1.0;
      index++;
    }
  }

  // use the permutation matrices to move the block to be sliced to the 
  // upper left corner of the matrix and return the block
  return (P1*A*P2).block(0,0,rowInds.size(), colInds.size());
}

// sorts a vector and returns the transformation
void getSortedIndices( const DoubleVector &x, UIntVector &res )
{
  // define the value and its original index as a vector of pairs
  vector<pair<Double,UInt> >dat(x.size());
  // assign the data
  for(UInt i=0; i<x.size(); i++) 
  {
    dat[i] = make_pair(x[i],i);
  }
  // sort according to the data
  sort(dat.begin(),dat.end());

  // get the result using the original indices in the pairs
  res = UIntVector(x.size());
  for(UInt i=0; i<x.size(); i++) 
  {
    res[i] = dat[i].second;
  }
}

// generates a julia file that can be used to generate a plot from physical DOF points and values 
void generateJulia2DRealPlot( vector<DoubleVector> pts, ComplexVector vals, string filename, string output_filename )
{
  // check if the input is consistent
  assert( pts.size() == vals.rows() );

  // open the .jl file that should be used
  ofstream file;
  file.open (filename);

  // include the PyPlot package
  file << "using PyPlot" << endl;

  // define the x-coordinates from the physical DOF points
  file << "x=[" << endl;
  for( UInt k=0; k<vals.rows(); k++ )
  {
    file << pts[k][0] << endl; 
  }
  file << "];" << endl;
  
  // define the y-coordinates from the physical DOF points
  file << "y=[" << endl;
  for( UInt k=0; k<vals.rows(); k++ )
  {
    file << pts[k][1] << endl; 
  }
  file << "];" << endl;

  // define the z-coordinates from thei real part of the values 
  file << "z=[" << endl;
  for( UInt k=0; k<vals.rows(); k++ )
  {
    file << real(vals[k]) << endl; 
  }
  file << "];" << endl;

  // the plotting command
  file << "tripcolor(x,y,z);" << endl;
  // enforce a colorbar
  file << "colorbar();" << endl;
  // use the given output filename in order to save the file in Julia as a .png
  file << "savefig(\"" << output_filename << "\");" << endl;
  // close the Julia file
  file.close();
}

