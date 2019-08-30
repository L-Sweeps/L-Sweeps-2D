#include "OutputVTK.h"
#include "Subdomain.h"
#include <iostream>
#include <fstream>

using namespace::std;

// The standard constructor
OutputVTK::OutputVTK(): Output()
{

}

// The copy constructor
OutputVTK::OutputVTK( const OutputVTK &other ): Output(other), mVals(other.mVals)
{

}

// The destructor
OutputVTK::~OutputVTK()
{

}

// The assignment operator
OutputVTK &OutputVTK::operator=( const OutputVTK &other )
{
  // only copy if other is different from this object
  if( this != &other )
  {
    // copy the object as an Output object
    (Output &)(*this) = (const Output&)other;

    // copy all member variables
    mVals = other.mVals;
  }
  return *this;
}

// The function that returns a pointer to a copy of the object
Output *OutputVTK::copy() const
{
  return new OutputVTK(*this);
}

// The function that generates the output file
void OutputVTK::print() const
{
  // determine the dimension
  UInt dim = mDD->getDim();

  // depending on the dimension, return the right auxiliary function
  // to write out the output file
  if( dim==2 )
    print2D();
  else
  {
    cerr << "ERROR: Only 2D meshes are supported to be printed in VTK format at tis point!" << endl;
    exit(1);
  }
}

// The function that sets the values
void OutputVTK::setVals( const vector<ComplexVector> &vals ) const
{
  // expand the values given to the functions
  mDD->expand( vals, mVals ); 
  // communicate the trace values so that 
  mDD->communicateTR2BLTraces( mVals );
  for( UInt i=0; i<mVals.size(); i++ )
    mDD->getCell(i)->expandDOFValsToPhysVals( mVals[i] );
}

vector<ComplexVector> OutputVTK::getVals() const
{
  return mVals;
}

// The function that writes a 2D output file
void OutputVTK::print2D() const
{
  // check if the dimension is correct
  assert( mDD->getDim()==2 );

  // allocate variables for the global points, values, and element info
  vector<DoubleVector> printPts(0);
  vector<Complex> printVals(0);
  vector<UIntVector> printEls(0);

  // get the number of total cells and loop over the cells to generate information
  UInt nCells = mDD->getNCells();
  for( UInt c=0; c<nCells; c++ )
  {
    // get the number of elements in the current cell
    UInt nEls = mDD->getCell(c)->getNEls();

    // allocate memory for a flag that determines which elements should be printed 
    vector<bool> temp(nEls);
    for( UInt i=0; i<nEls; i++ )
      temp[i] = false;
    
    // determine the elements that need to be printed 
    vector<UIntVector> elIndexSets;
    mDD->getCell(c)->getElIndexSets("LeftRight",false,elIndexSets);
    for( UInt i=3; i<=6; i++ )
    {
      for( UInt j=0; j<elIndexSets[i].size(); j++ )
      {
        temp[elIndexSets[i][j]] = true;   
      }
    }
    mDD->getCell(c)->getElIndexSets("UpDown",false,elIndexSets);
    for( UInt i=0; i<=2; i++ )
    {
      for( UInt j=0; j<elIndexSets[i].size(); j++ )
      {
        temp[elIndexSets[i][j]] = false;   
      }
    }
    for( UInt i=7; i<elIndexSets.size(); i++ )
    {
      for( UInt j=0; j<elIndexSets[i].size(); j++ )
      {
        temp[elIndexSets[i][j]] = false;   
      }
    }

    // generate a vector of element indices of the elements to be printed
    vector<UInt> elsToPrint(0);
    for( UInt i=0; i<temp.size(); i++ )
      if( temp[i] )
        elsToPrint.push_back(i);

    // get the local phsical points in the current cell
    vector<DoubleVector> pts;
    mDD->getCell(c)->getPhysPts( pts );

    // determine which points of pts are needed from this cell
    vector<bool> temp1(pts.size());
    for( UInt i=0; i<temp1.size(); i++ )
      temp1[i] = false;

    // expand printEls so that it includes the elements that need to be printed from this cell
    // for now from this cell, we store the local point indices for each element, they will be updated
    // to be global ones later
    UInt oldElPtSize = printEls.size();
    for( UInt i=0; i<elsToPrint.size(); i++ )
    {
      // get the local point indices
      UIntVector elPts;
      mDD->getCell(c)->getElPhysPtIndices( elsToPrint[i], elPts );
      // set the flags of the points that are need for the current element
      for( UInt j=0; j<elPts.size(); j++ )
        temp1[elPts[j]] = true;
      // add the element to printEls
      printEls.push_back(elPts);
    }

    // now add all values and points needed from this cell
    // this can be done using temp1 from above
    // As we add points, we also keep track of the indices that
    // map local point indices to global ones. This is stored in longToShortPts.
    vector<UInt> longToShortPts(temp1.size());
    for( UInt i=0; i<temp1.size(); i++ )
    {
      longToShortPts[i]=0;
      if( temp1[i] )
      {
        printVals.push_back( mVals[c][i] );
        printPts.push_back( pts[i] );
        longToShortPts[i] = printPts.size()-1;
      }
    }
    // now use longToShortPts to update the local point indices to global ones
    for( UInt i=oldElPtSize; i<printEls.size(); i++ )
      for( UInt j=0; j<printEls[i].size(); j++ )
        printEls[i][j] = longToShortPts[printEls[i][j]];
  }
  // all locally stored cells are processed now. We now have to communicate
  // data to get the global output file.

  // get the rank and number of processors
  UInt rank;
  MPI_Comm_rank( mDD->getComm(), (int*)&rank );
  UInt nprocs;
  MPI_Comm_size( mDD->getComm(), (int*)&nprocs );

  // if we are in rank 0, receive the data and add it to global variables
  // if we are in another rank, send the data
  if( rank==0 )
  {
    // loop over all ranks not equal to zero to receive the data
    for( UInt p=1; p<nprocs; p++ )
    {
      // get the size to be received from rank p
      MPI_Status status;
      MPI_Probe(p, 0, mDD->getComm(), &status);
      UInt size;
      MPI_Get_count( &status, MPI_DOUBLE, (int*)&size );

      // allocate memory for the x- and y-values of the points that need to be received
      DoubleVector printPtsXNew(size);
      DoubleVector printPtsYNew(size);
      // now receive the x- and y-values of the points
      MPI_Recv( printPtsXNew.data(), size, MPI_DOUBLE, p, 0, mDD->getComm(), &status );
      MPI_Recv( printPtsYNew.data(), size, MPI_DOUBLE, p, 1, mDD->getComm(), &status );
      // use the received information to add the new points
      UInt oldSize = printPts.size();
      printPts.resize( printPts.size()+size );
      for( UInt j=0; j<printPtsXNew.size(); j++ )
      {
        printPts[oldSize+j] = DoubleVector(2);
        printPts[oldSize+j][0] = printPtsXNew[j];
        printPts[oldSize+j][1] = printPtsYNew[j];
      }
      
      // get the number of elements that need to be received
      UInt numEls;
      MPI_Recv( &numEls, 1, MPI_INT, p, 0, mDD->getComm(), &status );
      // now add the received element information
      // allocate memory
      UInt oldNumEls = printEls.size();
      printEls.resize( printEls.size()+numEls );
      // loop over all elements, receive the information and update the information
      for( UInt j=0; j<numEls; j++ )
      {
        // get the number of point indices to be received for this element
        MPI_Probe(p, j, mDD->getComm(), &status);
        UInt size;
        MPI_Get_count( &status, MPI_INT, (int*)&size );
        // allocate memory for the point indices
        printEls[oldNumEls+j] = UIntVector(size);
        // receive the point indices
        MPI_Recv( printEls[oldNumEls+j].data(), size, MPI_INT, p, j, mDD->getComm(), &status );
        // update the received point indices to reflect the global point indices
        for( UInt i=0; i<printEls[oldNumEls+j].size(); i++ )
          printEls[oldNumEls+j][i] += oldSize;
      }

      // get the size of the values that need to be received from rank p
      MPI_Probe(p, 0, mDD->getComm(), &status);
      size = 0;
      MPI_Get_count( &status, MPI_DOUBLE, (int*)&size );

      // get the current size of the value vector
      oldSize = printVals.size();
      printVals.resize( printVals.size()+size);
      // receive the real and imaginary part of the values
      DoubleVector printValsNewRe(size);
      DoubleVector printValsNewIm(size);
      MPI_Recv( printValsNewRe.data(), size, MPI_DOUBLE, p, 0, mDD->getComm(), &status );
      MPI_Recv( printValsNewIm.data(), size, MPI_DOUBLE, p, 1, mDD->getComm(), &status );
      // use the received information to set the values   
      for( UInt i=0; i<size; i++ )
      {
        Complex im(0.0,1.0);
        printVals[oldSize+i] = printValsNewRe[i]+im*printValsNewIm[i];
      }
    }
  }
  else
  {
    // generate the informaiton that needs to be sent
    DoubleVector ptsXToSend( printPts.size() );
    DoubleVector ptsYToSend( printPts.size() );
    for( UInt i=0; i<printPts.size(); i++ )
    {
      ptsXToSend[i] = printPts[i][0];
      ptsYToSend[i] = printPts[i][1];
    }
    // send the information of the points
    MPI_Send( ptsXToSend.data(), ptsXToSend.size(), MPI_DOUBLE, 0, 0, mDD->getComm() );
    MPI_Send( ptsYToSend.data(), ptsYToSend.size(), MPI_DOUBLE, 0, 1, mDD->getComm() );

    // determine the number of elements that need to be sent
    UInt size = printEls.size();
    // send the number of elements
    MPI_Send( &size, 1, MPI_INT, 0, 0, mDD->getComm() );
    // loop over all elements and send the point index information
    for( UInt i=0; i<printEls.size(); i++ )
    {
      MPI_Send( printEls[i].data(), printEls[i].size(), MPI_INT, 0, i, mDD->getComm() );
    }

    // generate the informaiton that needs to be sent for each value (real and imaginary part)
    DoubleVector printValsRe( printVals.size() );
    DoubleVector printValsIm( printVals.size() );
    for( UInt i=0; i<printVals.size(); i++ )
    {
      printValsRe[i] = real( printVals[i] );
      printValsIm[i] = imag( printVals[i] );
    }
    // send the informaiton
    MPI_Send( printValsRe.data(), printValsRe.size(), MPI_DOUBLE, 0, 0, mDD->getComm() );
    MPI_Send( printValsIm.data(), printValsIm.size(), MPI_DOUBLE, 0, 1, mDD->getComm() );
  }

  // at this point, rank 0 holds all required global information, so we use it to write out the informaiton to the file
  if( rank==0 )
  {
    // number of total indices needed to describe the elements
    UInt totalElInds = 0;
    for( UInt i=0; i<printEls.size(); i++ )
      totalElInds += printEls[i].size()+1;

    // add the extension to the filename
    string filename = mFilename;
    filename.append(".vtk");
    // opent he file
    ofstream file;
    file.open(filename);

    // write the header for the output file
    file << "# vtk DataFile Version 1.0" << endl;
    file << "2D Unstuctured Grid of Linear Quads" << endl;
    file << "ASCII" << endl << endl;
    file << "DATASET UNSTRUCTURED_GRID" << endl;

    // write out all point info
    file << "POINTS " << printPts.size() << " float" << endl;
    for( UInt i=0; i<printPts.size(); i++ )
      file << printPts[i][0] << " " << printPts[i][1] << " 0.0" << endl;

    // write out all element info
    file << endl << "CELLS " << printEls.size() << " " << totalElInds << endl;
    for( UInt i=0; i<printEls.size(); i++ )
    {
      file << printEls[i].size();
      for( UInt j=0; j<printEls[i].size(); j++ )
        file << " " << printEls[i][j];
      file << endl;
    }
    
    // write out the types of elements
    file << endl << "CELL_TYPES " << printEls.size() << endl;
    for( UInt i=0; i<printEls.size(); i++ )
    {
      if( printEls[i].size()==3 )
        file << 5 << endl;
      else if( printEls[i].size()==4 )
        file << 9 << endl;
      else
      {
        cerr << "ERROR: Only trigs and quads are supported so far" << endl;
        exit(1);
      }
    }

    // write out the point data, we generate two variables,
    // one for the real and one for the imaginary part
    file << endl << "POINT_DATA " << printVals.size() << endl;
    file << "SCALARS wavefield_re float" << endl;
    file << "LOOKUP_TABLE default" << endl;
    Double val = 0.0;
    for( UInt i=0; i<printVals.size(); i++ )
      file << real(printVals[i]) << endl;
    file << "SCALARS wavefield_im float" << endl;
    file << "LOOKUP_TABLE default" << endl;
    for( UInt i=0; i<printVals.size(); i++ )
      file << imag(printVals[i]) << endl;

    // close the file
    file.close();
  }
}

