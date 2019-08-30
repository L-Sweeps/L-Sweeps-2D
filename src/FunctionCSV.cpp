#include "FunctionCSV.h"
#include <map>
#include <fstream>
#include <iostream>

using namespace std;

// The standard constructor, set up an empty object
FunctionCSV::FunctionCSV(): mA(0,0), mDomSet(false)
{

}

// The copy constructor
FunctionCSV::FunctionCSV( const FunctionCSV &other ): mFilename(other.mFilename), mFileDom(other.mFileDom), mFileN( other.mFileN), mH( other.mH ), mDom( other.mDom ), mDomA( other.mDomA ), mA( other.mA ), mDomSet(other.mDomSet)
{

}

// The destructor
FunctionCSV::~FunctionCSV()
{

}

// The assignment operator
FunctionCSV &FunctionCSV::operator=( const FunctionCSV &other )
{
  // Only copy if other is different from this object
  if( this != &other )
  {
    // copy all member variables
    mFilename = other.mFilename;
    mFileDom = other.mFileDom;
    mFileN = other.mFileN;
    mH = other.mH;
    mDom = other.mDom;
    mDomA = other.mDomA;
    mA = other.mA;
    mDomSet = other.mDomSet;
  }
  return *this;
}

// sets the filename and sets up the mesh corresponding to the file
void FunctionCSV::setFile( string filename, const vector<DoubleVector> &fileDom )
{
  // Set the filename
  mFilename = filename;

  // Set the domain corresponding to the file
  mFileDom = fileDom;


  // Now compute some initial information aboutthe mesh
  std::ifstream infile(filename);
  string line;

  // count the number of rows and columns in the matrix
  UInt nRows=0;
  UInt nCols=0;
  while( getline(infile,line) ) 
  {
    nRows++;

    istringstream ss(line);
    string val;

    nCols=0;
    while(std::getline(ss, val, ','))
      nCols++;
  }

  // Set the size of the mesh
  mFileN = UIntVector(2);
  mFileN[0] = nRows;
  mFileN[1] = nCols;

  // compute the mesh sizes from the size of the matrix
  mH = DoubleVector(2);
  mH[0] = (fileDom[0][1]-fileDom[0][0])/mFileN[1];
  mH[1] = (fileDom[1][1]-fileDom[1][0])/mFileN[0];
}

// This function sets the domain of definition
// The corresponding part of the matrix is read in from the file.
void FunctionCSV::setDomain( const vector<DoubleVector> &dom ) const
{
  // The domain of definition
  mDom = dom;
  // If the domain of definition is larger than the file corresponding to
  // the given file, truncate it. This allows for the correct extension of the
  // values in the normal direction
  mDom[0][0] = max( dom[0][0], mFileDom[0][0] );
  mDom[1][0] = max( dom[1][0], mFileDom[1][0] );
  mDom[0][1] = min( dom[0][1], mFileDom[0][1] );
  mDom[1][1] = min( dom[1][1], mFileDom[1][1] );

  // Find the first and last row of the mesh that is required to
  // cover the domain of definition
  UInt startRow = (mFileDom[1][1]-mDom[1][1])/mH[1];
  UInt endRow = (mFileDom[1][1]-mDom[1][0])/mH[1];
  // Find the first and last column of the mesh that is required to
  // cover the domain of definition
  UInt startCol = (mDom[0][0]-mFileDom[0][0])/mH[0];
  UInt endCol = (mDom[0][1]-mFileDom[0][0])/mH[0];
  endRow = min( endRow, mFileN[0]-1 );
  endCol = min( endCol, mFileN[1]-1 );

  // use the first/last row/column to define the domain of the
  // part of the mehs that covers the domain of definition
  mDomA = mDom;
  mDomA[0][0] = mFileDom[0][0]+startCol*mH[0];
  mDomA[0][1] = mFileDom[0][0]+(endCol+1)*mH[0];
  mDomA[1][0] = mFileDom[1][1]-(endRow+1)*mH[1];
  mDomA[1][1] = mFileDom[1][1]-startRow*mH[1];


  // now read in the part of the matrix
  // and determine the minimum/maximum values
  std::ifstream infile(mFilename);
  string line;
  UInt rowID = 0;

  // allocate memory for the matrix
  mA = ComplexMatrix(endRow-startRow+1,endCol-startCol+1);
  while( getline(infile,line) ) 
  {
    // if the row index is smaller than the starting row, just skip it
    if( rowID<startRow )
    {
      rowID++;
      continue;
    }
    // if the row index is larger than the end row, we are done and can break out of the loop
    if( rowID>endRow )
      break;

    // read in the line
    istringstream ss(line);
    string val;
    
    // loop over all elements in this line
    UInt colID = 0;
    while(std::getline(ss, val, ','))
    {
      // if the col index is smaller than the starting column, just skip it
      if( colID<startCol )
      {
        colID++;
        continue;
      }
      // if the col index is larger than the end col, we are done and can go to the next line
      if( colID>endCol )
        break;

      // convert the value of the wave speed to the squared slowness
      Complex tempVal = 1.0/pow(stof(val),2.0);
      
      // if it is the first value that is read in, set the minimum/maximum values accordingly
      // if not, then update the minimum/maximum value
      if(rowID==startRow && colID==startCol)
      {
        mMinReal = real(tempVal);
        mMaxReal = real(tempVal);
        mMinImag = imag(tempVal);
        mMaxImag = imag(tempVal);
      }
      else
      {
        mMinReal = min(mMinReal,real(tempVal));
        mMaxReal = max(mMaxReal,real(tempVal));
        mMinImag = min(mMinImag,imag(tempVal));
        mMaxImag = max(mMaxImag,imag(tempVal));
      }

      // set the value in the matrix
      mA(rowID-startRow,colID-startCol) = tempVal; 
      
      colID++;
    }
    rowID++;
  }

  // set the flag of a domain of definition to be set
  mDomSet=true;
}

// return the domain of definition
void FunctionCSV::getDomain( vector<DoubleVector> &res ) const
{
  res = mDom;
}

// the copy function
Function *FunctionCSV::copy() const
{
  return new FunctionCSV(*this);
}

// the function clearing the cache
void FunctionCSV::clearCache() const
{
  mDom.clear();
  mDomA.clear();
  mA = ComplexMatrix(0,0);
  mMinReal = 0.0;
  mMaxReal = 0.0;
  mMinImag = 0.0;
  mMaxImag = 0.0;
  mDomSet = false;
}

// the evaluation function
Complex FunctionCSV::eval( DoubleVector x ) const
{
  // make sure that the evaluation point has the same dimension as the domain
  assert( x.size()==mDom.size() );

  // if the point is outside the domain of definition, find the closest
  // point in the domain of definition in the normal direction
  for( UInt i=0; i<x.size(); i++ )
    if( x[i]<=mDomA[i][0]+0.5*mH[0] )
      x[i] = mDomA[i][0]+0.5*mH[0];
  for( UInt i=0; i<x.size(); i++ )
    if( x[i]>=mDomA[i][1]-0.5*mH[1] )
      x[i] = mDomA[i][1]-0.5*mH[1];

  // find the row and column index in the matrix
  // that corresponds to the given evaluation point
  UInt col = (x[0]-mDomA[0][0])/mH[0];
  UInt row = (mDomA[1][1]-x[1])/mH[1];

  // return the value in the matrix
  return mA(row,col);
}

// return the minimum real part of the value of the function 
Double FunctionCSV::getMinReal() const
{
  // if no domain of definition is set, we just use the entire domain of the file
  if( !mDomSet )
    setDomain(mFileDom);
  return mMinReal;  
}

// return the maximum real part of the value of the function 
Double FunctionCSV::getMaxReal() const
{
  // if no domain of definition is set, we just use the entire domain of the file
  if( !mDomSet )
    setDomain(mFileDom);
  return mMaxReal;  
}

// return the minimum imaginary part of the value of the function 
Double FunctionCSV::getMinImag() const
{
  // if no domain of definition is set, we just use the entire domain of the file
  if( !mDomSet )
    setDomain(mFileDom);
  return mMinImag;  
}

// return the maximum imaginary part of the value of the function 
Double FunctionCSV::getMaxImag() const
{
  // if no domain of definition is set, we just use the entire domain of the file
  if( !mDomSet )
    setDomain(mFileDom);
  return mMaxImag;  
}
