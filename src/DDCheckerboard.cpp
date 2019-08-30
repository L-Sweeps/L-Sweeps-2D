#include "DDCheckerboard.h"
#include "SubdomainStructured.h"
#include "LocalLinearSolverEigenSparseLU.h"
#include "utils.h"

#include <iostream>
#include <fstream>


using namespace::std;

// the copy constructor
DDCheckerboard::DDCheckerboard( const DDCheckerboard &other ): DD(other), mIndToProc(other.mIndToProc), mValsToKeep(other.mValsToKeep), mValsToKeepShifted(other.mValsToKeepShifted), mValsToKeepShiftedLR(other.mValsToKeepShiftedLR), mValsToKeepShiftedUD(other.mValsToKeepShiftedUD),mNCells1(other.mNCells1), mNCells2(other.mNCells2)
{
  // allocate memory for the cells
  mCells = vector<Subdomain*>(other.mCells.size());
  // loop over the cells and if there is a cell defined in other, copy it over
  for( UInt i=0; i<mCells.size(); i++ )
  {
    if( other.mCells[i]!=NULL )
      mCells[i] = other.mCells[i]->copy();
    else
      mCells[i] = NULL;
  }
}

// The copy function to properly copy the object
DD* DDCheckerboard::copy() const
{
  return new DDCheckerboard(*this);
}

// The destructor
DDCheckerboard::~DDCheckerboard()
{
  // loop over the cells and if they are defined, free the memory
  for( UInt i=0; i<mCells.size(); i++ )
  {
    if( mCells[i] != NULL )
      delete mCells[i];
  }
}

// the assignment operator
DDCheckerboard &DDCheckerboard::operator=( const DDCheckerboard &other )
{
  // only do something if other is not the same as this object
  if( this != &other )
  {
    // copy the communicator
    mComm = other.mComm;

    // allocate memory for the cells
    mCells = vector<Subdomain*>(other.mCells.size());
    // loop over the cells
    for( UInt i=0; i<mCells.size(); i++ )
    {
      // if a cell is defined, free the memory
      if( mCells[i]!=NULL )
        delete mCells[i];

      // if there is a cell defined at this cell, copy it over, otherwise set it to NULL
      if( other.mCells[i]!=NULL )
        mCells[i] = other.mCells[i]->copy();
      else
        mCells[i] = NULL;
    }
    
    // copy the processor assignment
    mIndToProc = other.mIndToProc;

    // copy the indices to keep
    mValsToKeep = other.mValsToKeep;
    mValsToKeepShifted = other.mValsToKeepShifted;
    mValsToKeepShiftedLR = other.mValsToKeepShiftedLR;
    mValsToKeepShiftedUD = other.mValsToKeepShiftedUD;

    // copy the number of cells
    mNCells1 = other.mNCells1;
    mNCells2 = other.mNCells2;
  }

  return *this;
}

// sets up a structured finite difference discretization
void DDCheckerboard::setUp(
    const Assembler* inputAssembler,
    const LocalLinearSolver* solver,
    vector<DoubleVector> domain,
    UIntVector nCells,
    UIntVector nElsInCell,
    MPI_Comm comm )
{
  // do some input checks
  assert( domain.size()==2 );
  assert( nCells.rows()==2 );

  // copy the assembler as a reference to set up the local problems
  Assembler *assembler = inputAssembler->copy();

  // get the order of discretization and determine the trace thickness
  UInt order = inputAssembler->getOrder();
  UInt traceThickness = int(ceil(order/2.0));

  // the number of elements inside one subdomain has to be sufficiently high in order
  // to accomodate all required (shifted) traces
  if( nElsInCell[0]<=4*traceThickness || nElsInCell[1]<=4*traceThickness )
  {
    cerr << "ERROR: Not enough elements in cell!" << endl;
    exit(1);
  }

  // assign the number of cells
  mNCells1 = nCells[0]; 
  mNCells2 = nCells[1]; 

  // number of elements in each direction (not including the PML)
  UInt n1 = nCells[0]*nElsInCell[0]+nCells[0]-1;
  UInt n2 = nCells[1]*nElsInCell[1]+nCells[1]-1;

  // the mesh size in each direction
  Double h1 = (domain[0][1]-domain[0][0])/n1;
  Double h2 = (domain[1][1]-domain[1][0])/n2;

  // we assign cells to processors in a row-based fashion
  // number of processors
  UInt nprocs;
  MPI_Comm_size(comm, (int*)&nprocs);

  // the local MPI rank
  UInt rank;
  MPI_Comm_rank(comm, (int*)&rank);

  // the minimum number of rows per processor
  UInt minRowsPerProc = UInt( floor(nCells[0]/nprocs) );

  // number of processors with more than the minimum amount of rows
  UInt numMaxRowProcs = nCells[0]-minRowsPerProc*nprocs;

  // number of processors with the minimum amount of rows
  UInt numMinRowProcs = nprocs-numMaxRowProcs;

  // allocate memory for the cells 
  if( rank<numMaxRowProcs )
  {
    mCells = vector<Subdomain*>(nCells[1]*(minRowsPerProc+1));
  }
  else
  {
    mCells = vector<Subdomain*>(nCells[1]*minRowsPerProc);
  }

  // compute the mapping from row to processor
  UIntVector rowToProc(nCells[0]);
  UInt index = 0;
  while( index<nCells[0] )
  {
    for( UInt r=0; r<nprocs; r++ )
    {
      if( index>=nCells[0] )
        break;
      rowToProc[index]=r;
      index++;
    }
  }
  MPI_Barrier( comm );

  // set up the mapping from subdomain to process
  mIndToProc = UIntMatrix(nCells[0],nCells[1]); 
  for( UInt i=0; i<nCells[0]; i++ )
  {
    for( UInt j=0; j<nCells[1]; j++ )
    {
      mIndToProc(i,j) = rowToProc[i];
    }
  }

  // allocate memory to store the indices of the local values
  mValsToKeep = vector<UIntVector>(mCells.size());
  mValsToKeepShifted = vector<UIntVector>(mCells.size());
  mValsToKeepShiftedLR = vector<UIntVector>(mCells.size());
  mValsToKeepShiftedUD = vector<UIntVector>(mCells.size());

  // loop over all cells and set up the subdomains
  for( UInt i=0; i<mCells.size()/nCells[1]; i++ )
  {
    for( UInt j=0; j<nCells[1]; j++ )
    {
      // assign the row index I
      UInt I=0;
      if( i<minRowsPerProc*nprocs )
        I=i*nprocs+rank;
      else
        I=numMinRowProcs*nprocs+rank;

      // assing the column index J
      UInt J=j;

      cout << "set up " << I << " / " << J << endl;

      // the computational domain (the domain including the PMLs) of the cell
      vector<DoubleVector> compDomain(2);
      compDomain[0] = DoubleVector(2);
      compDomain[1] = DoubleVector(2);

      // the limits of the traces
      DoubleVector locLims1(10);
      DoubleVector locLims2(10);

      // set up the computational domain and trace limits
      // the first direction
      // set the boundary of the subdomain on the left
      Double a = domain[0][0]+(nElsInCell[1]+1)*J*h1;

      // if we are in the first column, there are no traces on the left
      if( J>0 )
      {
      // Each limit should be in between two rows of degrees of freedom.
      // The limiters are defined to be left of the first row of DOFs
      // corresponding to
      // 0: trace0 non-shifted
      // 1: trace1 non-shifted
      // 2: trace0 shifted
      // 3: trace1 shifted
      // 4: remaining volume DOFs
        locLims1[0] = a-traceThickness*h1-0.5*h1; 
        locLims1[1] = a-0.5*h1; 
        locLims1[2] = a+traceThickness*h1-0.5*h1; 
        locLims1[3] = a+2*traceThickness*h1-0.5*h1; 
        locLims1[4] = a+3*traceThickness*h1-0.5*h1; 

        // extend the local domain to accomodate trace0 and keep a distance to the PML
        // to consistently capture information
        a -= 2*traceThickness*h1;
        if( a<domain[0][0] )
          a = domain[0][0];
      }
      // set the boundary of the subdomain on the right
      Double b = domain[0][0]+(nElsInCell[1]+1)*J*h1+nElsInCell[1]*h1;

      // if we are in the last column, there are no traces on the right 
      if( J<nCells[1]-1 )
      {
        // Each limit should be in between two rows of degrees of freedom.
        // The limiters are defined to be left of the first row of DOFs
        // corresponding to
        // 5: traceN non-shifted
        // 6: traceNP non-shifted
        // 7: traceN shifted
        // 8: traceNP shifted
        // 9: starting row of DOFs above traceNP shifted
        locLims1[5] = b-traceThickness*h1+0.5*h1;
        locLims1[6] = b+0.5*h1;
        locLims1[7] = b+traceThickness*h1+0.5*h1;
        locLims1[8] = b+2*traceThickness*h1+0.5*h1;
        locLims1[9] = b+3*traceThickness*h1+0.5*h1;

        // extend the local domain to accomodate traceNP, and the shifted traces
        // and keep a distance to the PML to consistently capture information
        b += 4*traceThickness*h1;
        if( b>domain[0][1] )
          b = domain[0][1];
      }
      // allocate memory for the domain of interest and set the boundaries accordingly
      vector<DoubleVector> intDomain(2);
      intDomain[0] = DoubleVector(2);
      intDomain[0][0] = a;
      intDomain[0][1] = b;


      // the second direction
 
      // set the boundary of the subdomain on the bottom
      a = domain[1][0]+(nElsInCell[0]+1)*I*h2;
      // if we are in the first row, there are no traces on the bottom
      if( I>0 )
      {
      // Each limit should be in between two rows of degrees of freedom.
      // The limiters are defined to be below the first row of DOFs
      // corresponding to
      // 0: trace0 non-shifted
      // 1: trace1 non-shifted
      // 2: trace0 shifted
      // 3: trace1 shifted
      // 4: remaining volume DOFs
        locLims2[0] = a-traceThickness*h2-0.5*h2;
        locLims2[1] = a-0.5*h2;
        locLims2[2] = a+traceThickness*h2-0.5*h2;
        locLims2[3] = a+2*traceThickness*h2-0.5*h2;
        locLims2[4] = a+3*traceThickness*h2-0.5*h2;
 
        // extend the local domain to accomodate trace0 and keep a distance to the PML
        // to consistently capture information
        a -= 2*traceThickness*h2;
        if( a<domain[1][0] )
          a = domain[1][0];
      }

      // set the boundary of the subdomain on the top
      b = domain[1][0]+(nElsInCell[0]+1)*I*h2+nElsInCell[0]*h2;

      // if we are in the last row, there are no traces on the right 
      if( I<nCells[0]-1 )
      {
        // Each limit should be in between two rows of degrees of freedom.
        // The limiters are defined to be below the first row of DOFs
        // corresponding to
        // 5: traceN non-shifted
        // 6: traceNP non-shifted
        // 7: traceN shifted
        // 8: traceNP shifted
        // 9: starting row of DOFs above traceNP shifted
        locLims2[5] = b-traceThickness*h2+0.5*h2;
        locLims2[6] = b+0.5*h2;
        locLims2[7] = b+traceThickness*h2+0.5*h2;
        locLims2[8] = b+2*traceThickness*h2+0.5*h2;
        locLims2[9] = b+3*traceThickness*h2+0.5*h2;
 
        // extend the local domain to accomodate traceNP, and the shifted traces
        // and keep a distance to the PML to consistently capture information
        b += 4*traceThickness*h2;
        if( b>domain[1][1] )
          b = domain[1][1];
      }

      // allocate memory for the domain of interest and set the boundaries accordingly
      intDomain[1] = DoubleVector(2);
      intDomain[1][0] = a;
      intDomain[1][1] = b;

      // determine the number of elements in each direction in the interior domain
      // of the local problem
      UInt n1 = (intDomain[0][1]-intDomain[0][0]+0.5*h1)/h1;
      UInt n2 = (intDomain[1][1]-intDomain[1][0]+0.5*h2)/h2;

      // set up the number of elements for input to the assembler
      UIntVector n(2);
      n[0] = n1;
      n[1] = n2;

      // define the subdomain
      SubdomainStructured subdomain;

      // assign the assembler to the subdomain
      subdomain.setAssembler(assembler);

      // set the interior domain in the subdomain
      subdomain.setDomain(intDomain);

      // set the number of elements
      // this is the number of elements in the interior, the PML region
      // is set up automatically
      subdomain.setN(n);

      // set the order of discretization
      subdomain.setOrder(order);

      // set the local solver
      subdomain.setLocalLinearSolver( solver );

      // set the row and column indices
      subdomain.setRowIndex( I );      
      subdomain.setColIndex( J );      

      // determine the minimum/maximum value in both directions
      // get the DOF points
      vector<DoubleVector> dofPts;
      subdomain.getDOFPts(dofPts);
      // use the DOF points to determin the minimum/maximum value
      Double minValX=0.0, maxValX=0.0;
      Double minValY=0.0, maxValY=0.0;
      for( UInt i=0; i<dofPts.size(); i++ )
      {
        if(dofPts[i][0]<minValX)
          minValX = dofPts[i][0];
        else if( dofPts[i][0]>maxValX )
          maxValX = dofPts[i][0];
        
        if(dofPts[i][1]<minValY)
          minValY = dofPts[i][1];
        else if( dofPts[i][1]>maxValY )
          maxValY = dofPts[i][1];
      }

      // use the minimum/maximum values to set up 
      // the limits at in the first/last row/column
      if( J==0 )
      {
        locLims1[0] = minValX-h1;
        locLims1[1] = minValX-h1;
        locLims1[2] = minValX-h1;
        locLims1[3] = minValX-h1;
        locLims1[4] = minValX-h1;
      }
      if( J==nCells[1]-1 )
      {
        locLims1[5] = maxValX+h1;
        locLims1[6] = maxValX+h1;
        locLims1[7] = maxValX+h1;
        locLims1[8] = maxValX+h1;
        locLims1[9] = maxValX+h1;
      }
      if( I==0 )
      {
        locLims2[0] = minValY-h2;
        locLims2[1] = minValY-h2;
        locLims2[2] = minValY-h2;
        locLims2[3] = minValY-h2;
        locLims2[4] = minValY-h2;
      }
      if( I==nCells[0]-1 )
      {
        locLims2[5] = maxValY+h2;
        locLims2[6] = maxValY+h2;
        locLims2[7] = maxValY+h2;
        locLims2[8] = maxValY+h2;
        locLims2[9] = maxValY+h2;
      }

      // set the limits for the definition of th traces in each direction
      vector<DoubleVector> lims(2);
      lims[0] = locLims1;
      lims[1] = locLims2;
      subdomain.setLimits( lims );

      // copy the subdomain into the cells
      mCells[i*nCells[1]+j] = subdomain.copy();

      // assemble the system matrix
      mCells[i*nCells[1]+j]->assembleSystemMatrix();


      // now find the indices that point to the values that should be stored
      // locally in this subdomain, this corresponds to all DOFs inside the 
      // domain of interest
      // 
      // get the indexSets in all directions and the physical DOF points
      vector<UIntVector> indexSets1, indexSets2;
      vector<UIntVector> indexSets1Shifted, indexSets2Shifted;
      mCells[i*nCells[1]+j]->getIndexSets("LeftRight",false, indexSets1 );
      mCells[i*nCells[1]+j]->getIndexSets("LeftRight",true, indexSets1Shifted );
      mCells[i*nCells[1]+j]->getIndexSets("UpDown",false, indexSets2 );
      mCells[i*nCells[1]+j]->getIndexSets("UpDown",true, indexSets2Shifted );
      vector<DoubleVector> pts;
      mCells[i*nCells[1]+j]->getDOFPts(pts);

      // define a bool vector as long as the number of DOFs in the subdomain
      // and denote by true or false if the DOF should be stored here or not
      // only trace1, interior volume DOFs, and traceN should be stored, so at first
      // use the LeftRight index sets to define those bool values
      // first for the non-shifted indices
      vector<bool> valsToKeep(subdomain.getSize());
      for(UInt k=0; k<indexSets1[0].size(); k++ )
        valsToKeep[indexSets1[0][k]] = false;
      for(UInt k=0; k<indexSets1[1].size(); k++ )
        valsToKeep[indexSets1[1][k]] = false;
      for(UInt k=0; k<indexSets1[2].size(); k++ )
        valsToKeep[indexSets1[2][k]] = true;
      for(UInt k=0; k<indexSets1[3].size(); k++ )
        valsToKeep[indexSets1[3][k]] = true;
      for(UInt k=0; k<indexSets1[4].size(); k++ )
        valsToKeep[indexSets1[4][k]] = true;
      for(UInt k=0; k<indexSets1[5].size(); k++ )
        valsToKeep[indexSets1[5][k]] = false;
      for(UInt k=0; k<indexSets1[6].size(); k++ )
        valsToKeep[indexSets1[6][k]] = false;

      vector<bool> valsToKeepShiftedUD(subdomain.getSize());
      for(UInt k=0; k<indexSets1[0].size(); k++ )
        valsToKeepShiftedUD[indexSets1[0][k]] = false;
      for(UInt k=0; k<indexSets1[1].size(); k++ )
        valsToKeepShiftedUD[indexSets1[1][k]] = false;
      for(UInt k=0; k<indexSets1[2].size(); k++ )
        valsToKeepShiftedUD[indexSets1[2][k]] = true;
      for(UInt k=0; k<indexSets1[3].size(); k++ )
        valsToKeepShiftedUD[indexSets1[3][k]] = true;
      for(UInt k=0; k<indexSets1[4].size(); k++ )
        valsToKeepShiftedUD[indexSets1[4][k]] = true;
      for(UInt k=0; k<indexSets1[5].size(); k++ )
        valsToKeepShiftedUD[indexSets1[5][k]] = false;
      for(UInt k=0; k<indexSets1[6].size(); k++ )
        valsToKeepShiftedUD[indexSets1[6][k]] = false;


      // then for the shifted indices
      vector<bool> valsToKeepShifted(subdomain.getSize());
      for(UInt k=0; k<indexSets1Shifted[0].size(); k++ )
        valsToKeepShifted[indexSets1Shifted[0][k]] = false;
      for(UInt k=0; k<indexSets1Shifted[1].size(); k++ )
        valsToKeepShifted[indexSets1Shifted[1][k]] = false;
      for(UInt k=0; k<indexSets1Shifted[2].size(); k++ )
        valsToKeepShifted[indexSets1Shifted[2][k]] = true;
      for(UInt k=0; k<indexSets1Shifted[3].size(); k++ )
        valsToKeepShifted[indexSets1Shifted[3][k]] = true;
      for(UInt k=0; k<indexSets1Shifted[4].size(); k++ )
        valsToKeepShifted[indexSets1Shifted[4][k]] = true;
      for(UInt k=0; k<indexSets1Shifted[5].size(); k++ )
        valsToKeepShifted[indexSets1Shifted[5][k]] = false;
      for(UInt k=0; k<indexSets1Shifted[6].size(); k++ )
        valsToKeepShifted[indexSets1Shifted[6][k]] = false;

      vector<bool> valsToKeepShiftedLR(subdomain.getSize());
      for(UInt k=0; k<indexSets1Shifted[0].size(); k++ )
        valsToKeepShiftedLR[indexSets1Shifted[0][k]] = false;
      for(UInt k=0; k<indexSets1Shifted[1].size(); k++ )
        valsToKeepShiftedLR[indexSets1Shifted[1][k]] = false;
      for(UInt k=0; k<indexSets1Shifted[2].size(); k++ )
        valsToKeepShiftedLR[indexSets1Shifted[2][k]] = true;
      for(UInt k=0; k<indexSets1Shifted[3].size(); k++ )
        valsToKeepShiftedLR[indexSets1Shifted[3][k]] = true;
      for(UInt k=0; k<indexSets1Shifted[4].size(); k++ )
        valsToKeepShiftedLR[indexSets1Shifted[4][k]] = true;
      for(UInt k=0; k<indexSets1Shifted[5].size(); k++ )
        valsToKeepShiftedLR[indexSets1Shifted[5][k]] = false;
      for(UInt k=0; k<indexSets1Shifted[6].size(); k++ )
        valsToKeepShiftedLR[indexSets1Shifted[6][k]] = false;



      // now we use the UpDown index sets to set the all index sets in this 
      // direction to false, except trace1, interior volume DOfs, and traceN
      // then valsToKeep is true exactly at the DOFs that should be stored here
      // first for the non-shifted indices
      for(UInt k=0; k<indexSets2[0].size(); k++ )
        valsToKeep[indexSets2[0][k]] = false;
      for(UInt k=0; k<indexSets2[1].size(); k++ )
        valsToKeep[indexSets2[1][k]] = false;
      for(UInt k=0; k<indexSets2[5].size(); k++ )
        valsToKeep[indexSets2[5][k]] = false;
      for(UInt k=0; k<indexSets2[6].size(); k++ )
        valsToKeep[indexSets2[6][k]] = false;

      for(UInt k=0; k<indexSets2[0].size(); k++ )
        valsToKeepShiftedLR[indexSets2[0][k]] = false;
      for(UInt k=0; k<indexSets2[1].size(); k++ )
        valsToKeepShiftedLR[indexSets2[1][k]] = false;
      for(UInt k=0; k<indexSets2[5].size(); k++ )
        valsToKeepShiftedLR[indexSets2[5][k]] = false;
      for(UInt k=0; k<indexSets2[6].size(); k++ )
        valsToKeepShiftedLR[indexSets2[6][k]] = false;


      // then for the shifted indices
      for(UInt k=0; k<indexSets2Shifted[0].size(); k++ )
        valsToKeepShifted[indexSets2Shifted[0][k]] = false;
      for(UInt k=0; k<indexSets2Shifted[1].size(); k++ )
        valsToKeepShifted[indexSets2Shifted[1][k]] = false;
      for(UInt k=0; k<indexSets2Shifted[5].size(); k++ )
        valsToKeepShifted[indexSets2Shifted[5][k]] = false;
      for(UInt k=0; k<indexSets2Shifted[6].size(); k++ )
        valsToKeepShifted[indexSets2Shifted[6][k]] = false;

      for(UInt k=0; k<indexSets2Shifted[0].size(); k++ )
        valsToKeepShiftedUD[indexSets2Shifted[0][k]] = false;
      for(UInt k=0; k<indexSets2Shifted[1].size(); k++ )
        valsToKeepShiftedUD[indexSets2Shifted[1][k]] = false;
      for(UInt k=0; k<indexSets2Shifted[5].size(); k++ )
        valsToKeepShiftedUD[indexSets2Shifted[5][k]] = false;
      for(UInt k=0; k<indexSets2Shifted[6].size(); k++ )
        valsToKeepShiftedUD[indexSets2Shifted[6][k]] = false;


      // count the number of DOFS that should be stored for the non-shifted indices
      UInt size=0;
      for( UInt k=0; k<valsToKeep.size(); k++ )
        if( valsToKeep[k] )
          size++;
     
      // allocate enough memory and set the indices that should be stored for the non-shifted indices
      mValsToKeep[i*nCells[1]+j] = UIntVector(size);
      UInt index = 0;
      for( UInt k=0; k<valsToKeep.size(); k++ )
      {
        if( valsToKeep[k] )
        {
          mValsToKeep[i*nCells[1]+j][index] = k;
          index++;
        }
      }

      // count the number of DOFS that should be stored for the shifted indices
      size=0;
      for( UInt k=0; k<valsToKeepShifted.size(); k++ )
        if( valsToKeepShifted[k] )
          size++;

      // allocate enough memory and set the indices that should be stored for the shifted indices
      mValsToKeepShifted[i*nCells[1]+j] = UIntVector(size);
      index = 0;
      for( UInt k=0; k<valsToKeepShifted.size(); k++ )
      {
        if( valsToKeepShifted[k] )
        {
          mValsToKeepShifted[i*nCells[1]+j][index] = k;
          index++;
        }
      }

      // count the number of DOFS that should be stored for the shifted indices
      size=0;
      for( UInt k=0; k<valsToKeepShiftedLR.size(); k++ )
        if( valsToKeepShiftedLR[k] )
          size++;

      // allocate enough memory and set the indices that should be stored for the shifted indices
      mValsToKeepShiftedLR[i*nCells[1]+j] = UIntVector(size);
      index = 0;
      for( UInt k=0; k<valsToKeepShiftedLR.size(); k++ )
      {
        if( valsToKeepShiftedLR[k] )
        {
          mValsToKeepShiftedLR[i*nCells[1]+j][index] = k;
          index++;
        }
      }

      // count the number of DOFS that should be stored for the shifted indices
      size=0;
      for( UInt k=0; k<valsToKeepShiftedUD.size(); k++ )
        if( valsToKeepShiftedUD[k] )
          size++;

      // allocate enough memory and set the indices that should be stored for the shifted indices
      mValsToKeepShiftedUD[i*nCells[1]+j] = UIntVector(size);
      index = 0;
      for( UInt k=0; k<valsToKeepShiftedUD.size(); k++ )
      {
        if( valsToKeepShiftedUD[k] )
        {
          mValsToKeepShiftedUD[i*nCells[1]+j][index] = k;
          index++;
        }
      }
    }
  }

  // delete the copy of the assembler used as a reference for setup
  delete assembler;

  // set the communicator
  mComm = comm;
  // wait until everything is done
  MPI_Barrier(comm);
}

// factorizes all subdomains
void DDCheckerboard::factorize()
{
  // get the rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);

  // loop over all subdomains and factorize them
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex();
    UInt J = mCells[i]->getColIndex();
    cout << "rank " << rank << ": (" << I << " , " << J << ")" << endl;
    mCells[i]->factorize();
  }
  
  // wait until everybody is done
  MPI_Barrier(mComm);
}

// assembles the rhs for a given RHS
void DDCheckerboard::assembleRHS( const RHS *rhs, vector<ComplexVector> &res ) const
{
  // allocate memory for the result
  vector<ComplexVector> resTemp(mCells.size());

  // loop over all local subdomains and assemble the entire rhs in the subdomain
  for(UInt i=0; i<mCells.size(); i++ )
  {
    mCells[i]->assembleRHS( rhs, resTemp[i] ); 
  }

  // reduce the local values to the ones that should be stored there
  reduce( resTemp, res );

  // wait until everybody is done
  MPI_Barrier(mComm);
}

// the function that applies the preconditioner (L-sweeps)
void DDCheckerboard::applyPrecond( const vector<ComplexVector> &x, vector<ComplexVector> &res ) const
{
  // separate the right hand side into the four windowed 
  // right hand sides.
  // The right-hand side for the non-shifted DD.
  vector<ComplexVector> rhs1;
  getSeparatedRHS(x,rhs1,false,false);
  // The right-hand side for the DD shifted in the Left-Right direction.
  vector<ComplexVector> rhs2;
  getSeparatedRHS(x,rhs2,true,false);
  // The right-hand side for the DD shifted in the Up-Down direction.
  vector<ComplexVector> rhs3;
  getSeparatedRHS(x,rhs3,false,true);
  // The right-hand side for the DD shifted in both directions.
  vector<ComplexVector> rhs4;
  getSeparatedRHS(x,rhs4,true,true);

  // apply the L-sweeps to the first DD and rhs
  vector<ComplexVector> res1,rhs;
  reduce( rhs1,rhs,false,false);
  applySweep( rhs,res1,false,false );
  
  // apply the L-sweeps to the second DD and rhs
  vector<ComplexVector> res2;
  reduce( rhs2,rhs,false,false);
  convertNonShiftedToShiftedLR( rhs, false );
  applySweep( rhs,res2,true,false );
  convertShiftedToNonShiftedLR( res2, false );
  
  // apply the L-sweeps to the third DD and rhs
  vector<ComplexVector> res3;
  reduce( rhs3,rhs,false,false);
  convertNonShiftedToShiftedUD( rhs, false );
  applySweep( rhs,res3,false,true );
  convertShiftedToNonShiftedUD( res3, false );
  
  // apply the L-sweeps to the fourth DD and rhs
  vector<ComplexVector> res4;
  reduce( rhs4,rhs,false,false);
  convertNonShiftedToShiftedLR( rhs, false );
  convertNonShiftedToShiftedUD( rhs, true );
  applySweep( rhs,res4,true,true );
  convertShiftedToNonShiftedUD( res4, true );
  convertShiftedToNonShiftedLR( res4, false );

  // add all results
  res = res1;
  for( UInt i=0; i<res.size(); i++ )
    res[i] += res2[i]+res3[i]+res4[i];

  // wait until everybody is done
  MPI_Barrier(mComm);
}

// the function that reduces local vectors that contain all DOFs in each local
// subdomain to vectors that only contain the DOFs that should be stored in each
// of the subdomains
void DDCheckerboard::reduce( const vector<ComplexVector> &x, vector<ComplexVector> &res, bool shiftedLR, bool shiftedUD ) const
{
  // allocate memory for the result
  res = vector<ComplexVector>(x.size());
  for( UInt i=0; i<res.size(); i++ )
  {
    // determine the right indices to keep
    UIntVector valsToKeep;
    if( shiftedLR && shiftedUD )
      valsToKeep = mValsToKeepShifted[i];
    else if( shiftedLR )
      valsToKeep = mValsToKeepShiftedLR[i];
    else if( shiftedUD )
      valsToKeep = mValsToKeepShiftedUD[i];
    else
      valsToKeep = mValsToKeep[i];

    // allocate the local memory
    res[i] = ComplexVector(valsToKeep.rows());

    // use the vector that holds the local indices that should be stored
    // in order to define the result
    for( UInt k=0; k<res[i].size(); k++ )
      res[i][k] = x[i][valsToKeep[k]];
  }
  
  // wait until everybody is done
  MPI_Barrier(mComm);
}

// the function that expands a reduced vector (one that only holds the values that
// should be stored in each subdomain) to a vector that holds all DOFs, the extension
// is by zero
void DDCheckerboard::expand( const vector<ComplexVector> &x, vector<ComplexVector> &res, bool shiftedLR, bool shiftedUD ) const
{
  // allocate the memory for the result
  res = vector<ComplexVector>(x.size());
  
  // loop over all subdomains
  for( UInt i=0; i<res.size(); i++ )
  {
    // determine the right indices to keep
    UIntVector valsToKeep;
    if( shiftedLR && shiftedUD )
      valsToKeep = mValsToKeepShifted[i];
    else if( shiftedLR )
      valsToKeep = mValsToKeepShiftedLR[i];
    else if( shiftedUD )
      valsToKeep = mValsToKeepShiftedUD[i];
    else
      valsToKeep = mValsToKeep[i];

    // allocate memory for the local subdomains
    res[i] = ComplexVector(mCells[i]->getSize());
    // set all values to zero
    for( UInt k=0; k<res[i].size(); k++ )
      res[i][k] = 0;
    // extract the values from the (reduced) input vector and put them at the
    // right locations
    for( UInt k=0; k<valsToKeep.size(); k++ )
      res[i][valsToKeep[k]] = x[i][k];
  }
  
  // wait until everybody is done
  MPI_Barrier(mComm);
}

// the function that applies the sweeps in all directions
void DDCheckerboard::applySweep( const vector<ComplexVector> &x, vector<ComplexVector> &res, bool shiftedLR, bool shiftedUD ) const
{
  // first expand the input vector
  vector<ComplexVector> xTemp;
  expand( x, xTemp, shiftedLR, shiftedUD );

  // allocate memory for the result and initialize it to zero
  vector<ComplexVector> resTemp(xTemp.size());
  for( UInt i=0; i<resTemp.size(); i++ )
    resTemp[i] = 0*xTemp[i];
  
  // define the vectors that will hold the trace values 
  vector<vector<ComplexVector>> tracesUD_right;
  vector<vector<ComplexVector>> tracesUD_left;
  vector<vector<ComplexVector>> tracesLR_up;
  vector<vector<ComplexVector>> tracesLR_down;
  vector<vector<ComplexVector>> tracesLR_loc;
  vector<vector<ComplexVector>> tracesUD_loc;

  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);
  if( rank==0 )
  cout << "stage 1" << endl;
  // add the contributions of the local solutions and extract the local traces
  addLocSol( xTemp, resTemp, tracesLR_loc, tracesUD_loc, shiftedLR, shiftedUD );
  
  if( rank==0 )
  cout << "stage 2" << endl;
  if( rank==0 )
  cout << "right" << endl;
  // add the contributions of the right sweep and extract the local traces in the Up-Down direction
  addRightSweep( resTemp, tracesLR_loc, tracesUD_right, shiftedLR, shiftedUD );
  if( rank==0 )
  cout << "left" << endl;
  // add the contributions of the left sweep and extract the local traces in the Up-Down direction
  addLeftSweep( resTemp, tracesLR_loc, tracesUD_left, shiftedLR, shiftedUD );
  if( rank==0 )
  cout << "up" << endl;
  // add the contributions of the up sweep and extract the local traces in the Left-Right direction
  addUpSweep( resTemp, tracesLR_up, tracesUD_loc, shiftedLR, shiftedUD );
  if( rank==0 )
  cout << "down" << endl;
  // add the contributions of the down sweep and extract the local traces in the Left-Right direction
  addDownSweep( resTemp, tracesLR_down, tracesUD_loc, shiftedLR, shiftedUD );
 
  if( rank==0 )
  cout << "stage 3" << endl;
  if( rank==0 )
  cout << "BL2TR" << endl;
  // add the contributions from the sweep from botton left to top right
  addBL2TRSweep( resTemp, tracesLR_up, tracesUD_right, shiftedLR, shiftedUD );
  if( rank==0 )
  cout << "TR2BL" << endl;
  // add the contributions from the sweep from top right to bottom left
  addTR2BLSweep( resTemp, tracesLR_down, tracesUD_left, shiftedLR, shiftedUD );
  if( rank==0 )
  cout << "BR2TL" << endl;
  // add the contributions from the sweep from bottom right to top left
  addBR2TLSweep( resTemp, tracesLR_up, tracesUD_left, shiftedLR, shiftedUD );
  if( rank==0 )
  cout << "TL2BR" << endl;
  // add the contributions from the sweep from top left to bottom right
  addTL2BRSweep( resTemp, tracesLR_down, tracesUD_right, shiftedLR, shiftedUD );

  // for now we return the result as it is (in unreduced form)
  // but later the result should be reduced to only the DOFs that should be
  // stored locally
  reduce( resTemp, res, shiftedLR, shiftedUD );
}

// the function that adds the constributions of the local solutions (stage 1) 
void DDCheckerboard::addLocSol( const vector<ComplexVector> &x, vector<ComplexVector> &res, vector<vector<ComplexVector>> &tracesLR, vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);

  // allocate memory for the UpDown traces (those are the ones that will
  // be returned)
  tracesUD = vector<vector<ComplexVector>>(4);
  tracesUD[0] = vector<ComplexVector>(x.size());
  tracesUD[1] = vector<ComplexVector>(x.size());
  tracesUD[2] = vector<ComplexVector>(x.size());
  tracesUD[3] = vector<ComplexVector>(x.size());

  // allocate memory for the LeftRight traces (those are just used in this
  // function and will not be returned)
  tracesLR = vector<vector<ComplexVector>>(4);
  tracesLR[0] = vector<ComplexVector>(x.size());
  tracesLR[1] = vector<ComplexVector>(x.size());
  tracesLR[2] = vector<ComplexVector>(x.size());
  tracesLR[3] = vector<ComplexVector>(x.size());

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2;

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<mCells.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = d;
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();

    // get the index sets for the local subdomain
    vector<UIntVector> indexSets1, indexSets2;
    mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
    mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);

    // now solve for the contributions of the local sources
    //
    // define the temporary result
    ComplexVector resTemp, rhsTemp(x[index].size());

    // define the rhs
    for( UInt k=0; k<rhsTemp.size(); k++ )
      rhsTemp[k]=x[index][k];

    // solve the local problem
    mCells[index]->solve( rhsTemp, resTemp );

    // add the contribution of the local sources to the result
    res[index] += resTemp;

    // use the solution to define the local trace DOFs that are required
    //
    // allocate memory for the UpDown traces and set the to zero
    // (no information for the UpSweep has to be stored here becuase this
    // will be taken care of by the up/down Sweep)
    tracesUD[0][index] = ComplexVector(indexSets2[1].rows());
    tracesUD[1][index] = ComplexVector(indexSets2[2].rows());
    tracesUD[2][index] = ComplexVector(indexSets2[4].rows());
    tracesUD[3][index] = ComplexVector(indexSets2[5].rows());
    for( UInt k=0; k<indexSets2[1].rows(); k++ )
      tracesUD[0][index][k] = resTemp[indexSets2[1][k]]; 
    for( UInt k=0; k<indexSets2[2].rows(); k++ )
      tracesUD[1][index][k] = resTemp[indexSets2[2][k]]; 
    for( UInt k=0; k<indexSets2[4].rows(); k++ )
      tracesUD[2][index][k] = resTemp[indexSets2[4][k]]; 
    for( UInt k=0; k<indexSets2[5].rows(); k++ )
      tracesUD[3][index][k] = resTemp[indexSets2[5][k]]; 
    // allocate memory for the LeftRight traces and set them using the solution
    // note, we only set the traces to the right because these are just used in this
    // function and we only compute the right sweep here
    tracesLR[0][index] = ComplexVector(indexSets1[1].rows());
    tracesLR[1][index] = ComplexVector(indexSets1[2].rows());
    tracesLR[2][index] = ComplexVector(indexSets1[4].rows());
    tracesLR[3][index] = ComplexVector(indexSets1[5].rows());
    for( UInt k=0; k<indexSets1[1].rows(); k++ )
      tracesLR[0][index][k] = resTemp[indexSets1[1][k]];
    for( UInt k=0; k<indexSets1[2].rows(); k++ )
      tracesLR[1][index][k] = resTemp[indexSets1[2][k]];
    for( UInt k=0; k<indexSets1[4].rows(); k++ )
      tracesLR[2][index][k] = resTemp[indexSets1[4][k]];
    for( UInt k=0; k<indexSets1[5].rows(); k++ )
      tracesLR[3][index][k] = resTemp[indexSets1[5][k]];
  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}

// the function that adds the left sweep (part of stage 2)
void DDCheckerboard::addLeftSweep( vector<ComplexVector> &res, const vector<vector<ComplexVector>> &tracesLR, vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);

  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  // note that we define the diagonals from top right to bottom left
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    I = mNCells1-1-I;
    J = mNCells2-1-J;

    diagInds[i] = I+J;
    
    I = mNCells1-1-I;
    J = mNCells2-1-J;
    
    locInds(I,J)=i;
  }
 
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 
  
  // allocate memory for the UpDown traces (those are the ones that will
  // be returned)
  tracesUD = vector<vector<ComplexVector>>(4);
  tracesUD[0] = vector<ComplexVector>(mCells.size());
  tracesUD[1] = vector<ComplexVector>(mCells.size());
  tracesUD[2] = vector<ComplexVector>(mCells.size());
  tracesUD[3] = vector<ComplexVector>(mCells.size());

  // allocate memory for the LeftRight traces (those are just used in this
  // function and will not be returned)
  vector<vector<ComplexVector>> tracesLRTemp = tracesLR;

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2;

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
    
    // get the index sets for the local subdomain
    vector<UIntVector> indexSets1, indexSets2;
    mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
    mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);

    // use the solution to define the local trace DOFs that are required
    //
    // allocate memory for the UpDown traces and set the to zero
    // (no information for the UpSweep has to be stored here becuase this
    // will be taken care of by the up/down Sweep)
    tracesUD[0][index] = ComplexVector(indexSets2[1].rows());
    tracesUD[1][index] = ComplexVector(indexSets2[2].rows());
    tracesUD[2][index] = ComplexVector(indexSets2[4].rows());
    tracesUD[3][index] = ComplexVector(indexSets2[5].rows());
    for( UInt k=0; k<indexSets2[1].rows(); k++ )
      tracesUD[0][index][k] = 0; 
    for( UInt k=0; k<indexSets2[2].rows(); k++ )
      tracesUD[1][index][k] = 0; 
    for( UInt k=0; k<indexSets2[4].rows(); k++ )
      tracesUD[2][index][k] = 0; 
    for( UInt k=0; k<indexSets2[5].rows(); k++ )
      tracesUD[3][index][k] = 0; 

    // if we are not in the last column, we also have to add the contribution
    // that comes from the information coming into the domain from the element
    // on the right
    if( JJ<mNCells2-1 )
    {
      // get the trace information from the right element
      //
      // if the right element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ+1) )
      {
        UInt index_p1 = locInds(II,JJ+1);
        tracesLRTemp[2][index] = tracesLRTemp[0][index_p1];
        tracesLRTemp[3][index] = tracesLRTemp[1][index_p1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesLR2_re(indexSets1[4].size());   
        DoubleVector tracesLR2_im(indexSets1[4].size());   
        DoubleVector tracesLR3_re(indexSets1[5].size());   
        DoubleVector tracesLR3_im(indexSets1[5].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesLR2_re.data(), indexSets1[4].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesLR2_im.data(), indexSets1[4].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesLR3_re.data(), indexSets1[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesLR3_im.data(), indexSets1[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesLRTemp[2][index] = tracesLR2_re+im*tracesLR2_im;
        tracesLRTemp[3][index] = tracesLR3_re+im*tracesLR3_im;
      }
     
      // use the local system matrix and the trace information to define the
      // equivalent sources on the left
      vector<ComplexVector> equivSrc;
      mCells[index]->getEquivalentSources("Right",shiftedLR,shiftedUD,tracesLRTemp[3][index],tracesLRTemp[2][index],equivSrc);

      // now we can solve the local system using the equivalent sources
      // set the rhs and allocate memory for the result
      ComplexVector rhsTemp(mCells[index]->getSize());
      for( UInt i=0; i<rhsTemp.size(); i++ )
        rhsTemp[i] = 0;
      ComplexVector resTemp = rhsTemp;
      for( UInt k=0; k<indexSets1[4].rows(); k++ )
        rhsTemp[indexSets1[4][k]] = equivSrc[1][k];
      for( UInt k=0; k<indexSets1[5].rows(); k++ )
        rhsTemp[indexSets1[5][k]] = equivSrc[0][k];

      // solve the system
      mCells[index]->solve( rhsTemp, resTemp );
      // add the contribtuion
      res[index] += resTemp;

     // use the solution to define the local trace DOFs that are required
     // we do need to store the UpDown trace information from the left sweep because
     // this will not be handled by the up/down sweep
     // again, we only store the left traces for the LeftRight traces becuase
     // we are only considering the right sweep here
      for( UInt k=0; k<indexSets2[1].rows(); k++ )
        tracesUD[0][index][k] += resTemp[indexSets2[1][k]];
      for( UInt k=0; k<indexSets2[2].rows(); k++ )
        tracesUD[1][index][k] += resTemp[indexSets2[2][k]];
      for( UInt k=0; k<indexSets2[4].rows(); k++ )
        tracesUD[2][index][k] += resTemp[indexSets2[4][k]];
      for( UInt k=0; k<indexSets2[5].rows(); k++ )
        tracesUD[3][index][k] += resTemp[indexSets2[5][k]];
      for( UInt k=0; k<indexSets1[1].rows(); k++ )
        tracesLRTemp[0][index][k] += resTemp[indexSets1[1][k]];
      for( UInt k=0; k<indexSets1[2].rows(); k++ )
        tracesLRTemp[1][index][k] += resTemp[indexSets1[2][k]];
    }

    // if we are not in the first column and the element on the left
    // is not in the same process, we need to send the trace info to the left
    if( JJ>0 && mIndToProc(II,JJ)!=mIndToProc(II,JJ-1) )
    {
      // devide up the traceN info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR0_re(tracesLRTemp[0][index].size());
      DoubleVector tracesLR0_im(tracesLRTemp[0][index].size());
      for( UInt k=0; k<tracesLR0_re.size(); k++ )
      {
        tracesLR0_re[k] = real(tracesLRTemp[0][index][k]);
        tracesLR0_im[k] = imag(tracesLRTemp[0][index][k]);
      }
      // send the info
      MPI_Send( tracesLR0_re.data(), tracesLR0_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+0, mComm );
      MPI_Send( tracesLR0_im.data(), tracesLR0_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+1, mComm );

      // devide up the traceNP info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR1_re(tracesLRTemp[1][index].size());
      DoubleVector tracesLR1_im(tracesLRTemp[1][index].size());
      for( UInt k=0; k<tracesLR1_re.size(); k++ )
      {
        tracesLR1_re[k] = real(tracesLRTemp[1][index][k]);
        tracesLR1_im[k] = imag(tracesLRTemp[1][index][k]);
      }
      // send the info
      MPI_Send( tracesLR1_re.data(), tracesLR1_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+2, mComm );
      MPI_Send( tracesLR1_im.data(), tracesLR1_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+3, mComm );
    }
  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}


// the function that adds the constribution of the right sweep (part of stage 2)
void DDCheckerboard::addRightSweep( vector<ComplexVector> &res, const vector<vector<ComplexVector>> &tracesLR, vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);

  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    diagInds[i] = I+J;
    locInds(I,J)=i;
  }
  
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 
  
  // allocate memory for the UpDown traces (those are the ones that will
  // be returned)
  tracesUD = vector<vector<ComplexVector>>(4);
  tracesUD[0] = vector<ComplexVector>(mCells.size());
  tracesUD[1] = vector<ComplexVector>(mCells.size());
  tracesUD[2] = vector<ComplexVector>(mCells.size());
  tracesUD[3] = vector<ComplexVector>(mCells.size());

  vector<vector<ComplexVector>> tracesLRTemp = tracesLR;

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2;

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
    
    // get the index sets for the local subdomain
    vector<UIntVector> indexSets1, indexSets2;
    mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
    mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);
   
    // use the solution to define the local trace DOFs that are required
    //
    // allocate memory for the UpDown traces and set the to zero
    // (no information for the UpSweep has to be stored here becuase this
    // will be taken care of by the up/down Sweep)
    tracesUD[0][index] = ComplexVector(indexSets2[1].rows());
    tracesUD[1][index] = ComplexVector(indexSets2[2].rows());
    tracesUD[2][index] = ComplexVector(indexSets2[4].rows());
    tracesUD[3][index] = ComplexVector(indexSets2[5].rows());
    for( UInt k=0; k<indexSets2[1].rows(); k++ )
      tracesUD[0][index][k] = 0; 
    for( UInt k=0; k<indexSets2[2].rows(); k++ )
      tracesUD[1][index][k] = 0; 
    for( UInt k=0; k<indexSets2[4].rows(); k++ )
      tracesUD[2][index][k] = 0; 
    for( UInt k=0; k<indexSets2[5].rows(); k++ )
      tracesUD[3][index][k] = 0; 

    // if we are not in the first column, we also have to add the contribution
    // that comes from the information coming into the domain from the element
    // on the left
    if( JJ>0 )
    {
      // get the trace information from the left element
      //
      // if the left element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ-1) )
      {
        UInt index_m1 = locInds(II,JJ-1);
        tracesLRTemp[0][index] = tracesLRTemp[2][index_m1];
        tracesLRTemp[1][index] = tracesLRTemp[3][index_m1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesLR0_re(indexSets1[1].size());   
        DoubleVector tracesLR0_im(indexSets1[1].size());   
        DoubleVector tracesLR1_re(indexSets1[2].size());   
        DoubleVector tracesLR1_im(indexSets1[2].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesLR0_re.data(), indexSets1[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesLR0_im.data(), indexSets1[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesLR1_re.data(), indexSets1[2].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesLR1_im.data(), indexSets1[2].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesLRTemp[0][index] = tracesLR0_re+im*tracesLR0_im;
        tracesLRTemp[1][index] = tracesLR1_re+im*tracesLR1_im;
      }

      // use the local system matrix and the trace information to define the
      // equivalent sources on the left
      vector<ComplexVector> equivSrc;
      mCells[index]->getEquivalentSources("Left",shiftedLR,shiftedUD,tracesLRTemp[0][index],tracesLRTemp[1][index],equivSrc);

      // now we can solve the local system using the equivalent sources
      // set the rhs and allocate memory for the result
      ComplexVector rhsTemp(mCells[index]->getSize());
      for( UInt i=0; i<rhsTemp.size(); i++ )
        rhsTemp[i] = 0;
      ComplexVector resTemp = rhsTemp;
      assert( equivSrc.size()==2 );
      assert( equivSrc[0].size()==indexSets1[1].size() );
      assert( equivSrc[1].size()==indexSets1[2].size() );
      for( UInt k=0; k<indexSets1[1].rows(); k++ )
        rhsTemp[indexSets1[1][k]] = equivSrc[0][k];
      for( UInt k=0; k<indexSets1[2].rows(); k++ )
        rhsTemp[indexSets1[2][k]] = equivSrc[1][k];
      
      // solve the system
      mCells[index]->solve( rhsTemp, resTemp );
      // add the contribution
      res[index] += resTemp;

     // use the solution to define the local trace DOFs that are required
     // we do need to store the UpDown trace information from the right sweep because
     // this will not be handled by the up sweep
     // again, we only store the right traces for the LeftRight traces becuase
     // we are only considering the right sweep here
      for( UInt k=0; k<indexSets2[1].rows(); k++ )
        tracesUD[0][index][k] += resTemp[indexSets2[1][k]];
      for( UInt k=0; k<indexSets2[2].rows(); k++ )
        tracesUD[1][index][k] += resTemp[indexSets2[2][k]];
      for( UInt k=0; k<indexSets2[4].rows(); k++ )
        tracesUD[2][index][k] += resTemp[indexSets2[4][k]];
      for( UInt k=0; k<indexSets2[5].rows(); k++ )
        tracesUD[3][index][k] += resTemp[indexSets2[5][k]];
      for( UInt k=0; k<indexSets1[4].rows(); k++ )
        tracesLRTemp[2][index][k] += resTemp[indexSets1[4][k]];
      for( UInt k=0; k<indexSets1[5].rows(); k++ )
        tracesLRTemp[3][index][k] += resTemp[indexSets1[5][k]];
    }

    // if we are not in the last column and the element on the right
    // is not in the same process, we need to send the trace info to the right
    if( JJ<mNCells2-1 && mIndToProc(II,JJ)!=mIndToProc(II,JJ+1) )
    {
      // devide up the traceN info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR2_re(tracesLRTemp[2][index].size());
      DoubleVector tracesLR2_im(tracesLRTemp[2][index].size());
      for( UInt k=0; k<tracesLR2_re.size(); k++ )
      {
        tracesLR2_re[k] = real(tracesLRTemp[2][index][k]);
        tracesLR2_im[k] = imag(tracesLRTemp[2][index][k]);
      }
      // send the info
      MPI_Send( tracesLR2_re.data(), tracesLR2_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+0, mComm );
      MPI_Send( tracesLR2_im.data(), tracesLR2_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+1, mComm );

      // devide up the traceNP info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR3_re(tracesLRTemp[3][index].size());
      DoubleVector tracesLR3_im(tracesLRTemp[3][index].size());
      for( UInt k=0; k<tracesLR3_re.size(); k++ )
      {
        tracesLR3_re[k] = real(tracesLRTemp[3][index][k]);
        tracesLR3_im[k] = imag(tracesLRTemp[3][index][k]);
      }
      // send the info
      MPI_Send( tracesLR3_re.data(), tracesLR3_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+2, mComm );
      MPI_Send( tracesLR3_im.data(), tracesLR3_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+3, mComm );
    }

  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}

// the function that adds the contribution of the down sweep (stage 2)
void DDCheckerboard::addDownSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> &tracesLR, const vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);
 
  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    I = mNCells1-1-I;
    J = mNCells2-1-J;
    diagInds[i] = I+J;
    
    I = mNCells1-1-I;
    J = mNCells2-1-J;
    locInds(I,J)=i;
  }
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 

  // allocate memory for the UpDown traces (those are just used in this
  // function and will not be returned)
  vector<vector<ComplexVector>> tracesUDTemp = tracesUD;
  
  // allocate memory for the LeftRight traces (those are the ones that will
  // be returned)
  tracesLR = vector<vector<ComplexVector>>(4);
  tracesLR[0] = vector<ComplexVector>(mCells.size());
  tracesLR[1] = vector<ComplexVector>(mCells.size());
  tracesLR[2] = vector<ComplexVector>(mCells.size());
  tracesLR[3] = vector<ComplexVector>(mCells.size());

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2; 

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
  
    // get the index sets for the local subdomain
    vector<UIntVector> indexSets1, indexSets2;
    mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
    mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);

    // use the solution to define the local trace DOFs that are required
    //
    // allocate memory for the LeftRgiht traces and set the to zero
    // (no information for the LeftRight has to be stored here becuase this
    // will be taken care of by the right/left Sweep)
    tracesLR[0][index] = ComplexVector(indexSets1[1].rows());
    tracesLR[1][index] = ComplexVector(indexSets1[2].rows());
    tracesLR[2][index] = ComplexVector(indexSets1[4].rows());
    tracesLR[3][index] = ComplexVector(indexSets1[5].rows());
    for( UInt k=0; k<indexSets1[1].rows(); k++ )
      tracesLR[0][index][k] = 0; 
    for( UInt k=0; k<indexSets1[2].rows(); k++ )
      tracesLR[1][index][k] = 0; 
    for( UInt k=0; k<indexSets1[4].rows(); k++ )
      tracesLR[2][index][k] = 0; 
    for( UInt k=0; k<indexSets1[5].rows(); k++ )
      tracesLR[3][index][k] = 0; 
    
    // if we are not in the last row, we also have to add the contribution
    // that comes from the information coming into the domain from the element
    // below
    if( II<mNCells1-1 )
    {
      // get the trace information from the top element
      //
      // if the top element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II+1,JJ) )
      {
        UInt index_p1 = locInds(II+1,JJ);
        tracesUDTemp[2][index] = tracesUDTemp[0][index_p1];
        tracesUDTemp[3][index] = tracesUDTemp[1][index_p1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesUD2_re(indexSets2[4].size());   
        DoubleVector tracesUD2_im(indexSets2[4].size());   
        DoubleVector tracesUD3_re(indexSets2[5].size());   
        DoubleVector tracesUD3_im(indexSets2[5].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesUD2_re.data(), indexSets2[4].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesUD2_im.data(), indexSets2[4].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesUD3_re.data(), indexSets2[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesUD3_im.data(), indexSets2[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesUDTemp[2][index] = tracesUD2_re+im*tracesUD2_im;
        tracesUDTemp[3][index] = tracesUD3_re+im*tracesUD3_im;
      }

      // use the local system matrix and the trace information to define the
      // equivalent sources on the bottom 
      vector<ComplexVector> equivSrc;
      mCells[index]->getEquivalentSources("Top",shiftedLR,shiftedUD,tracesUDTemp[3][index],tracesUDTemp[2][index],equivSrc);

      // now we can solve the local system using the equivalent sources
      // set the rhs and allocate memory for the result
      ComplexVector rhsTemp(mCells[index]->getSize());
      for( UInt i=0; i<rhsTemp.size(); i++ )
        rhsTemp[i] = 0;
      ComplexVector resTemp = rhsTemp;
      assert( equivSrc.size()==2 );
      assert( equivSrc[0].size()==indexSets2[5].size() );
      assert( equivSrc[1].size()==indexSets2[4].size() );
      for( UInt k=0; k<indexSets2[4].rows(); k++ )
        rhsTemp[indexSets2[4][k]] = equivSrc[1][k];
      for( UInt k=0; k<indexSets2[5].rows(); k++ )
        rhsTemp[indexSets2[5][k]] = equivSrc[0][k];

      // solve the system
      mCells[index]->solve( rhsTemp, resTemp );
      // add the contribution
      res[index] += resTemp;
      
     // use the solution to define the local trace DOFs that are required
     // we do need to store the LeftRight trace information from the down sweep because
     // this will not be handled by the right/left sweep
     // again, we only store the bottom traces for the UpDown traces becuase
     // we are only considering the up sweep here
      for( UInt k=0; k<indexSets1[1].rows(); k++ )
        tracesLR[0][index][k] += resTemp[indexSets1[1][k]];
      for( UInt k=0; k<indexSets1[2].rows(); k++ )
        tracesLR[1][index][k] += resTemp[indexSets1[2][k]];
      for( UInt k=0; k<indexSets1[4].rows(); k++ )
        tracesLR[2][index][k] += resTemp[indexSets1[4][k]];
      for( UInt k=0; k<indexSets1[5].rows(); k++ )
        tracesLR[3][index][k] += resTemp[indexSets1[5][k]];
      for( UInt k=0; k<indexSets2[1].rows(); k++ )
        tracesUDTemp[0][index][k] += resTemp[indexSets2[1][k]];
      for( UInt k=0; k<indexSets2[2].rows(); k++ )
        tracesUDTemp[1][index][k] += resTemp[indexSets2[2][k]];
      }

    // if we are not in the first row and the element on the bottom
    // is not in the same process, we need to send the trace info to the bottom
    if( II>0 && mIndToProc(II,JJ)!=mIndToProc(II-1,JJ) )
    {
      // devide up the trace0 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD0_re(tracesUDTemp[0][index].size());
      DoubleVector tracesUD0_im(tracesUDTemp[0][index].size());
      for( UInt k=0; k<tracesUD0_re.size(); k++ )
      {
        tracesUD0_re[k] = real(tracesUDTemp[0][index][k]);
        tracesUD0_im[k] = imag(tracesUDTemp[0][index][k]);
      }
      // send the info
      MPI_Send( tracesUD0_re.data(), tracesUD0_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+0, mComm );
      MPI_Send( tracesUD0_im.data(), tracesUD0_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+1, mComm );

      // devide up the trace1 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD1_re(tracesUDTemp[1][index].size());
      DoubleVector tracesUD1_im(tracesUDTemp[1][index].size());
      for( UInt k=0; k<tracesUD1_re.size(); k++ )
      {
        tracesUD1_re[k] = real(tracesUDTemp[1][index][k]);
        tracesUD1_im[k] = imag(tracesUDTemp[1][index][k]);
      }
      // send the info
      MPI_Send( tracesUD1_re.data(), tracesUD1_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+2, mComm );
      MPI_Send( tracesUD1_im.data(), tracesUD1_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+3, mComm );
    }
  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}


// the function that adds the constribution of the up sweep (part of stage 2)
void DDCheckerboard::addUpSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> &tracesLR, const vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);
 
  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    diagInds[i] = I+J;
    locInds(I,J)=i;
  }
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 

  vector<vector<ComplexVector>> tracesUDTemp = tracesUD;
  
  // allocate memory for the LeftRight traces (those are the ones that will
  // be returned)
  tracesLR = vector<vector<ComplexVector>>(4);
  tracesLR[0] = vector<ComplexVector>(mCells.size());
  tracesLR[1] = vector<ComplexVector>(mCells.size());
  tracesLR[2] = vector<ComplexVector>(mCells.size());
  tracesLR[3] = vector<ComplexVector>(mCells.size());

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2; 

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
   
    // get the index sets for the local subdomain
    vector<UIntVector> indexSets1, indexSets2;
    mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
    mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);

    // use the solution to define the local trace DOFs that are required
    //
    // allocate memory for the LeftRgiht traces and set the to zero
    // (no information for the LeftRight has to be stored here becuase this
    // will be taken care of by the right/left Sweep)
    tracesLR[0][index] = ComplexVector(indexSets1[1].rows());
    tracesLR[1][index] = ComplexVector(indexSets1[2].rows());
    tracesLR[2][index] = ComplexVector(indexSets1[4].rows());
    tracesLR[3][index] = ComplexVector(indexSets1[5].rows());
    for( UInt k=0; k<indexSets1[1].rows(); k++ )
      tracesLR[0][index][k] = 0; 
    for( UInt k=0; k<indexSets1[2].rows(); k++ )
      tracesLR[1][index][k] = 0; 
    for( UInt k=0; k<indexSets1[4].rows(); k++ )
      tracesLR[2][index][k] = 0; 
    for( UInt k=0; k<indexSets1[5].rows(); k++ )
      tracesLR[3][index][k] = 0; 
    
    // if we are not in the first row, we also have to add the contribution
    // that comes from the information coming into the domain from the element
    // below
    if( II>0 )
    {
      // get the trace information from the bottom element
      //
      // if the bottom element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II-1,JJ) )
      {
        UInt index_m1 = locInds(II-1,JJ);
        tracesUDTemp[0][index] = tracesUDTemp[2][index_m1];
        tracesUDTemp[1][index] = tracesUDTemp[3][index_m1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesUD0_re(indexSets2[1].size());   
        DoubleVector tracesUD0_im(indexSets2[1].size());   
        DoubleVector tracesUD1_re(indexSets2[2].size());   
        DoubleVector tracesUD1_im(indexSets2[2].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesUD0_re.data(), indexSets2[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesUD0_im.data(), indexSets2[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesUD1_re.data(), indexSets2[2].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesUD1_im.data(), indexSets2[2].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesUDTemp[0][index] = tracesUD0_re+im*tracesUD0_im;
        tracesUDTemp[1][index] = tracesUD1_re+im*tracesUD1_im;
      }

      // use the local system matrix and the trace information to define the
      // equivalent sources on the bottom 
      vector<ComplexVector> equivSrc;
      mCells[index]->getEquivalentSources("Bottom",shiftedLR,shiftedUD,tracesUDTemp[0][index],tracesUDTemp[1][index],equivSrc);

      // now we can solve the local system using the equivalent sources
      // set the rhs and allocate memory for the result
      ComplexVector rhsTemp(mCells[index]->getSize());
      for( UInt i=0; i<rhsTemp.size(); i++ )
        rhsTemp[i] = 0;
      ComplexVector resTemp = rhsTemp;
      assert( equivSrc.size()==2 );
      assert( equivSrc[0].size()==indexSets2[1].size() );
      assert( equivSrc[1].size()==indexSets2[2].size() );
      for( UInt k=0; k<indexSets2[1].rows(); k++ )
        rhsTemp[indexSets2[1][k]] = equivSrc[0][k];
      for( UInt k=0; k<indexSets2[2].rows(); k++ )
        rhsTemp[indexSets2[2][k]] = equivSrc[1][k];
     
      // solve the system
      mCells[index]->solve( rhsTemp, resTemp );
      // add the contribution
      res[index] += resTemp;
      
     // use the solution to define the local trace DOFs that are required
     // we do need to store the LeftRight trace information from the up sweep because
     // this will not be handled by the right sweep
     // again, we only store the top traces for the UpDown traces becuase
     // we are only considering the up sweep here
      for( UInt k=0; k<indexSets1[1].rows(); k++ )
        tracesLR[0][index][k] += resTemp[indexSets1[1][k]];
      for( UInt k=0; k<indexSets1[2].rows(); k++ )
        tracesLR[1][index][k] += resTemp[indexSets1[2][k]];
      for( UInt k=0; k<indexSets1[4].rows(); k++ )
        tracesLR[2][index][k] += resTemp[indexSets1[4][k]];
      for( UInt k=0; k<indexSets1[5].rows(); k++ )
        tracesLR[3][index][k] += resTemp[indexSets1[5][k]];
      for( UInt k=0; k<indexSets2[4].rows(); k++ )
        tracesUDTemp[2][index][k] += resTemp[indexSets2[4][k]];
      for( UInt k=0; k<indexSets2[5].rows(); k++ )
        tracesUDTemp[3][index][k] += resTemp[indexSets2[5][k]];
      }

    // if we are not in the last row and the element on the top
    // is not in the same process, we need to send the trace info to the top
    if( II<mNCells1-1 && mIndToProc(II,JJ)!=mIndToProc(II+1,JJ) )
    {
      // devide up the traceN info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD2_re(tracesUDTemp[2][index].size());
      DoubleVector tracesUD2_im(tracesUDTemp[2][index].size());
      for( UInt k=0; k<tracesUD2_re.size(); k++ )
      {
        tracesUD2_re[k] = real(tracesUDTemp[2][index][k]);
        tracesUD2_im[k] = imag(tracesUDTemp[2][index][k]);
      }
      // send the info
      MPI_Send( tracesUD2_re.data(), tracesUD2_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+0, mComm );
      MPI_Send( tracesUD2_im.data(), tracesUD2_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+1, mComm );

      // devide up the traceNP info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD3_re(tracesUDTemp[3][index].size());
      DoubleVector tracesUD3_im(tracesUDTemp[3][index].size());
      for( UInt k=0; k<tracesUD3_re.size(); k++ )
      {
        tracesUD3_re[k] = real(tracesUDTemp[3][index][k]);
        tracesUD3_im[k] = imag(tracesUDTemp[3][index][k]);
      }
      // send the info
      MPI_Send( tracesUD3_re.data(), tracesUD3_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+2, mComm );
      MPI_Send( tracesUD3_im.data(), tracesUD3_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+3, mComm );
    }
  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}

// the functiont that adds the contribution of the diagonal sweep from top left to bottom right (part of stage 3)
void DDCheckerboard::addTL2BRSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> tracesLR, vector<vector<ComplexVector>> tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);
  
  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    I = mNCells1-1-I;
    diagInds[i] = I+J;

    I = mNCells1-1-I;
    locInds(I,J)=i;
  }
 
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2; 

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
   
    // only do something if we are not in the last row or the first column
    // (these cases were handled by the right/left and up/down sweep)
    if( II<mNCells1-1 && JJ>0 )
    {
      // get the index sets for the local subdomain
      vector<UIntVector> indexSets1, indexSets2;
      mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
      mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);

      // get the trace information from the left and top element
      //
      // if the left element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ-1) )
      {
        UInt index_m1 = locInds(II,JJ-1);
        tracesLR[0][index] = tracesLR[2][index_m1];
        tracesLR[1][index] = tracesLR[3][index_m1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesLR0_re(indexSets1[1].size());   
        DoubleVector tracesLR0_im(indexSets1[1].size());   
        DoubleVector tracesLR1_re(indexSets1[2].size());   
        DoubleVector tracesLR1_im(indexSets1[2].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesLR0_re.data(), indexSets1[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesLR0_im.data(), indexSets1[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesLR1_re.data(), indexSets1[2].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesLR1_im.data(), indexSets1[2].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesLR[0][index] = tracesLR0_re+im*tracesLR0_im;
        tracesLR[1][index] = tracesLR1_re+im*tracesLR1_im;
      }
      // if the top element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II+1,JJ) )
      {
        UInt index_p1 = locInds(II+1,JJ);
        tracesUD[2][index] = tracesUD[0][index_p1];
        tracesUD[3][index] = tracesUD[1][index_p1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesUD2_re(indexSets2[4].size());   
        DoubleVector tracesUD2_im(indexSets2[4].size());   
        DoubleVector tracesUD3_re(indexSets2[5].size());   
        DoubleVector tracesUD3_im(indexSets2[5].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesUD2_re.data(), indexSets2[4].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesUD2_im.data(), indexSets2[4].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesUD3_re.data(), indexSets2[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesUD3_im.data(), indexSets2[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesUD[2][index] = tracesUD2_re+im*tracesUD2_im;
        tracesUD[3][index] = tracesUD3_re+im*tracesUD3_im;
      }
  
      // extract the L-shaped traces
      vector<ComplexVector> LTraceVals;
      mCells[index]->extractLTraces( "TopLeft", shiftedLR, shiftedUD, tracesLR[0][index], tracesLR[1][index], tracesUD[3][index], tracesUD[2][index], LTraceVals );

      // compute the equivalent sources on the boundary (the right-hand side contributions to compute the polarized wave field)
      vector<ComplexVector> equivSrcs;
      mCells[index]->getEquivalentSources( "TopLeft", shiftedLR, shiftedUD, LTraceVals[0], LTraceVals[1], equivSrcs );

      // get the indicex sets corresponding to L-shaped traces
      vector<UIntVector> LTraceIndexSets;
      mCells[index]->getLTraceIndexSets( "TopLeft", shiftedLR, shiftedUD, LTraceIndexSets );

      // now we can solve the local system using the equivalent sources
      // set the rhs and allocate memory for the result
      ComplexVector rhsTemp(mCells[index]->getSize());
      for( UInt i=0; i<rhsTemp.size(); i++ )
        rhsTemp[i] = 0;
      ComplexVector resTemp = rhsTemp;
      for( UInt k=0; k<LTraceIndexSets[0].rows(); k++ )
        rhsTemp[LTraceIndexSets[0][k]] = equivSrcs[0][k];
      for( UInt k=0; k<LTraceIndexSets[1].rows(); k++ )
        rhsTemp[LTraceIndexSets[1][k]] = equivSrcs[1][k];

      // solve the system
      mCells[index]->solve( rhsTemp, resTemp );
      // add the sweep contribution to the result
      res[index] += resTemp;
      // add the traces that are not interior DOFs to the current solution
      // so that the traces are approriately extracted
      for( UInt k=0; k<indexSets1[1].size(); k++ )
        resTemp[indexSets1[1][k]] = tracesLR[0][index][k];
      for( UInt k=0; k<indexSets2[5].size(); k++ )
        resTemp[indexSets2[5][k]] = tracesUD[3][index][k];

     // use the solution to define the local trace DOFs that are required
     // we do need to store both the UpDown and LeftRight trace information
     // because both are needed for the L-sweep.
     // However, we only need to store the right and bottom traces because
     // those are the only ones needed in the next steps of this L sweep
      for( UInt k=0; k<indexSets1[4].rows(); k++ )
        tracesLR[2][index][k] += resTemp[indexSets1[4][k]];
      for( UInt k=0; k<indexSets1[5].rows(); k++ )
        tracesLR[3][index][k] += resTemp[indexSets1[5][k]];
      for( UInt k=0; k<indexSets2[1].rows(); k++ )
        tracesUD[0][index][k] += resTemp[indexSets2[1][k]];
      for( UInt k=0; k<indexSets2[2].rows(); k++ )
        tracesUD[1][index][k] += resTemp[indexSets2[2][k]];
    }

    // if we are not in the last column or last row and the element on the right
    // is not in the same process, we need to send the trace info to the right
    if( JJ<mNCells2-1 && mIndToProc(II,JJ)!=mIndToProc(II,JJ+1) && II<mNCells1-1 )
    {
      // devide up the trace0 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR2_re(tracesLR[2][index].size());
      DoubleVector tracesLR2_im(tracesLR[2][index].size());
      for( UInt k=0; k<tracesLR2_re.size(); k++ )
      {
        tracesLR2_re[k] = real(tracesLR[2][index][k]);
        tracesLR2_im[k] = imag(tracesLR[2][index][k]);
      }
      // send the info
      MPI_Send( tracesLR2_re.data(), tracesLR2_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+0, mComm );
      MPI_Send( tracesLR2_im.data(), tracesLR2_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+1, mComm );

      // devide up the trace1 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR3_re(tracesLR[3][index].size());
      DoubleVector tracesLR3_im(tracesLR[3][index].size());
      for( UInt k=0; k<tracesLR3_re.size(); k++ )
      {
        tracesLR3_re[k] = real(tracesLR[3][index][k]);
        tracesLR3_im[k] = imag(tracesLR[3][index][k]);
      }
      // send the info
      MPI_Send( tracesLR3_re.data(), tracesLR3_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+2, mComm );
      MPI_Send( tracesLR3_im.data(), tracesLR3_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+3, mComm );
    }

    // if we are not in the first row or first column and the element below
    // is not in the same process, we need to send the trace info to the bottom
    if( II>0 && mIndToProc(II,JJ)!=mIndToProc(II-1,JJ) && JJ>0 )
    {
      // devide up the trace0 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD0_re(tracesUD[0][index].size());
      DoubleVector tracesUD0_im(tracesUD[0][index].size());
      for( UInt k=0; k<tracesUD0_re.size(); k++ )
      {
        tracesUD0_re[k] = real(tracesUD[0][index][k]);
        tracesUD0_im[k] = imag(tracesUD[0][index][k]);
      }
      // send the info
      MPI_Send( tracesUD0_re.data(), tracesUD0_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+0, mComm );
      MPI_Send( tracesUD0_im.data(), tracesUD0_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+1, mComm );

      // devide up the trace1 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD1_re(tracesUD[1][index].size());
      DoubleVector tracesUD1_im(tracesUD[1][index].size());
      for( UInt k=0; k<tracesUD1_re.size(); k++ )
      {
        tracesUD1_re[k] = real(tracesUD[1][index][k]);
        tracesUD1_im[k] = imag(tracesUD[1][index][k]);
      }
      // send the info
      MPI_Send( tracesUD1_re.data(), tracesUD1_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+2, mComm );
      MPI_Send( tracesUD1_im.data(), tracesUD1_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+3, mComm );
    }
  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}

// the function that adds the contribution of the diagonal sweep from top right to bottom left
void DDCheckerboard::addTR2BLSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> tracesLR, vector<vector<ComplexVector>> tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);
  
  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    I = mNCells1-1-I;
    J = mNCells2-1-J;
    diagInds[i] = I+J;

    I = mNCells1-1-I;
    J = mNCells2-1-J;
    locInds(I,J)=i;
  }
 
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2; 

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
   
    // only do something if we are not in the last row or the last column
    // (these cases were handled by the right/left and up/down sweep)
    if( II<mNCells1-1 && JJ<mNCells2-1 )
    {
      // get the index sets for the local subdomain
      vector<UIntVector> indexSets1, indexSets2;
      mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
      mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);

      // get the trace information from the right and top element
      //
      // if the right element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ+1) )
      {
        UInt index_p1 = locInds(II,JJ+1);
        tracesLR[2][index] = tracesLR[0][index_p1];
        tracesLR[3][index] = tracesLR[1][index_p1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesLR2_re(indexSets1[4].size());   
        DoubleVector tracesLR2_im(indexSets1[4].size());   
        DoubleVector tracesLR3_re(indexSets1[5].size());   
        DoubleVector tracesLR3_im(indexSets1[5].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesLR2_re.data(), indexSets1[4].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesLR2_im.data(), indexSets1[4].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesLR3_re.data(), indexSets1[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesLR3_im.data(), indexSets1[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesLR[2][index] = tracesLR2_re+im*tracesLR2_im;
        tracesLR[3][index] = tracesLR3_re+im*tracesLR3_im;
      }
      // if the top element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II+1,JJ) )
      {
        UInt index_p1 = locInds(II+1,JJ);
        tracesUD[2][index] = tracesUD[0][index_p1];
        tracesUD[3][index] = tracesUD[1][index_p1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesUD2_re(indexSets2[4].size());   
        DoubleVector tracesUD2_im(indexSets2[4].size());   
        DoubleVector tracesUD3_re(indexSets2[5].size());   
        DoubleVector tracesUD3_im(indexSets2[5].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesUD2_re.data(), indexSets2[4].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesUD2_im.data(), indexSets2[4].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesUD3_re.data(), indexSets2[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesUD3_im.data(), indexSets2[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesUD[2][index] = tracesUD2_re+im*tracesUD2_im;
        tracesUD[3][index] = tracesUD3_re+im*tracesUD3_im;
      }
      
      // extract the L-shaped traces
      vector<ComplexVector> LTraceVals;
      mCells[index]->extractLTraces( "TopRight", shiftedLR, shiftedUD, tracesLR[3][index], tracesLR[2][index], tracesUD[3][index], tracesUD[2][index], LTraceVals );

      // compute the equivalent sources on the boundary (the right-hand side contributions to compute the polarized wave field)
      vector<ComplexVector> equivSrcs;
      mCells[index]->getEquivalentSources( "TopRight", shiftedLR, shiftedUD, LTraceVals[0], LTraceVals[1], equivSrcs );

      // get the indicex sets corresponding to L-shaped traces
      vector<UIntVector> LTraceIndexSets;
      mCells[index]->getLTraceIndexSets( "TopRight", shiftedLR, shiftedUD, LTraceIndexSets );

      // now we can solve the local system using the equivalent sources
      // set the rhs and allocate memory for the result
      ComplexVector rhsTemp(mCells[index]->getSize());
      for( UInt i=0; i<rhsTemp.size(); i++ )
        rhsTemp[i] = 0;
      ComplexVector resTemp = rhsTemp;
      for( UInt k=0; k<LTraceIndexSets[0].rows(); k++ )
        rhsTemp[LTraceIndexSets[0][k]] = equivSrcs[0][k];
      for( UInt k=0; k<LTraceIndexSets[1].rows(); k++ )
        rhsTemp[LTraceIndexSets[1][k]] = equivSrcs[1][k];

      // solve the system
      mCells[index]->solve( rhsTemp, resTemp );
      // add the sweep contribution to the result
      res[index] += resTemp;
      // add the traces that are not interior DOFs to the current solution
      // so that the traces are approriately extracted
      for( UInt k=0; k<indexSets1[5].size(); k++ )
        resTemp[indexSets1[5][k]] = tracesLR[3][index][k];
      for( UInt k=0; k<indexSets2[5].size(); k++ )
        resTemp[indexSets2[5][k]] = tracesUD[3][index][k];

     // use the solution to define the local trace DOFs that are required
     // we do need to store both the UpDown and LeftRight trace information
     // because both are needed for the L-sweep.
     // However, we only need to store the left and bottom traces because
     // those are the only ones needed in the next steps of this L sweep
      for( UInt k=0; k<indexSets1[1].rows(); k++ )
        tracesLR[0][index][k] += resTemp[indexSets1[1][k]];
      for( UInt k=0; k<indexSets1[2].rows(); k++ )
        tracesLR[1][index][k] += resTemp[indexSets1[2][k]];
      for( UInt k=0; k<indexSets2[1].rows(); k++ )
        tracesUD[0][index][k] += resTemp[indexSets2[1][k]];
      for( UInt k=0; k<indexSets2[2].rows(); k++ )
        tracesUD[1][index][k] += resTemp[indexSets2[2][k]];
    }

    // if we are not in the first column or last row and the element on the left
    // is not in the same process, we need to send the trace info to the left
    if( JJ>0 && mIndToProc(II,JJ)!=mIndToProc(II,JJ-1) && II<mNCells1-1 )
    {
      // devide up the trace0 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR0_re(tracesLR[0][index].size());
      DoubleVector tracesLR0_im(tracesLR[0][index].size());
      for( UInt k=0; k<tracesLR0_re.size(); k++ )
      {
        tracesLR0_re[k] = real(tracesLR[0][index][k]);
        tracesLR0_im[k] = imag(tracesLR[0][index][k]);
      }
      // send the info
      MPI_Send( tracesLR0_re.data(), tracesLR0_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+0, mComm );
      MPI_Send( tracesLR0_im.data(), tracesLR0_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+1, mComm );

      // devide up the trace1 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR1_re(tracesLR[1][index].size());
      DoubleVector tracesLR1_im(tracesLR[1][index].size());
      for( UInt k=0; k<tracesLR1_re.size(); k++ )
      {
        tracesLR1_re[k] = real(tracesLR[1][index][k]);
        tracesLR1_im[k] = imag(tracesLR[1][index][k]);
      }
      // send the info
      MPI_Send( tracesLR1_re.data(), tracesLR1_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+2, mComm );
      MPI_Send( tracesLR1_im.data(), tracesLR1_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+3, mComm );
    }

    // if we are not in the first row or last column and the element below
    // is not in the same process, we need to send the trace info to the bottom
    if( II>0 && mIndToProc(II,JJ)!=mIndToProc(II-1,JJ) && JJ<mNCells2-1 )
    {
      // devide up the trace0 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD0_re(tracesUD[0][index].size());
      DoubleVector tracesUD0_im(tracesUD[0][index].size());
      for( UInt k=0; k<tracesUD0_re.size(); k++ )
      {
        tracesUD0_re[k] = real(tracesUD[0][index][k]);
        tracesUD0_im[k] = imag(tracesUD[0][index][k]);
      }
      // send the info
      MPI_Send( tracesUD0_re.data(), tracesUD0_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+0, mComm );
      MPI_Send( tracesUD0_im.data(), tracesUD0_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+1, mComm );

      // devide up the trace1 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD1_re(tracesUD[1][index].size());
      DoubleVector tracesUD1_im(tracesUD[1][index].size());
      for( UInt k=0; k<tracesUD1_re.size(); k++ )
      {
        tracesUD1_re[k] = real(tracesUD[1][index][k]);
        tracesUD1_im[k] = imag(tracesUD[1][index][k]);
      }
      // send the info
      MPI_Send( tracesUD1_re.data(), tracesUD1_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+2, mComm );
      MPI_Send( tracesUD1_im.data(), tracesUD1_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+3, mComm );
    }
  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}

// the functiont that adds the contribution of the diagonal sweep from bottom right to top left
void DDCheckerboard::addBR2TLSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> tracesLR, vector<vector<ComplexVector>> tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);
  
  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    J = mNCells2-1-J;
    diagInds[i] = I+J;

    J = mNCells2-1-J;
    locInds(I,J)=i;
  }
 
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2; 

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
   
    // only do something if we are not in the first row or the last column
    // (these cases were handled by the right/left and up/down sweep)
    if( II>0 && JJ<mNCells2-1 )
    {
      // get the index sets for the local subdomain
      vector<UIntVector> indexSets1, indexSets2;
      mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
      mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);

      // get the trace information from the right and bottom element
      //
      // if the right element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ+1) )
      {
        UInt index_p1 = locInds(II,JJ+1);
        tracesLR[2][index] = tracesLR[0][index_p1];
        tracesLR[3][index] = tracesLR[1][index_p1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesLR2_re(indexSets1[4].size());   
        DoubleVector tracesLR2_im(indexSets1[4].size());   
        DoubleVector tracesLR3_re(indexSets1[5].size());   
        DoubleVector tracesLR3_im(indexSets1[5].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesLR2_re.data(), indexSets1[4].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesLR2_im.data(), indexSets1[4].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesLR3_re.data(), indexSets1[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesLR3_im.data(), indexSets1[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesLR[2][index] = tracesLR2_re+im*tracesLR2_im;
        tracesLR[3][index] = tracesLR3_re+im*tracesLR3_im;
      }
      // if the bottom element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II-1,JJ) )
      {
        UInt index_m1 = locInds(II-1,JJ);
        tracesUD[0][index] = tracesUD[2][index_m1];
        tracesUD[1][index] = tracesUD[3][index_m1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesUD0_re(indexSets2[1].size());   
        DoubleVector tracesUD0_im(indexSets2[1].size());   
        DoubleVector tracesUD1_re(indexSets2[2].size());   
        DoubleVector tracesUD1_im(indexSets2[2].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesUD0_re.data(), indexSets2[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesUD0_im.data(), indexSets2[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesUD1_re.data(), indexSets2[2].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesUD1_im.data(), indexSets2[2].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesUD[0][index] = tracesUD0_re+im*tracesUD0_im;
        tracesUD[1][index] = tracesUD1_re+im*tracesUD1_im;
      }

      // extract the L-shaped traces
      vector<ComplexVector> LTraceVals;
      mCells[index]->extractLTraces( "BottomRight", shiftedLR, shiftedUD, tracesLR[3][index], tracesLR[2][index], tracesUD[0][index], tracesUD[1][index], LTraceVals );

      // compute the equivalent sources on the boundary (the right-hand side contributions to compute the polarized wave field)
      vector<ComplexVector> equivSrcs;
      mCells[index]->getEquivalentSources( "BottomRight", shiftedLR, shiftedUD, LTraceVals[0], LTraceVals[1], equivSrcs );

      // get the indicex sets corresponding to L-shaped traces
      vector<UIntVector> LTraceIndexSets;
      mCells[index]->getLTraceIndexSets( "BottomRight", shiftedLR, shiftedUD, LTraceIndexSets );

      // now we can solve the local system using the equivalent sources
      // set the rhs and allocate memory for the result
      ComplexVector rhsTemp(mCells[index]->getSize());
      for( UInt i=0; i<rhsTemp.size(); i++ )
        rhsTemp[i] = 0;
      ComplexVector resTemp = rhsTemp;
      for( UInt k=0; k<LTraceIndexSets[0].rows(); k++ )
        rhsTemp[LTraceIndexSets[0][k]] = equivSrcs[0][k];
      for( UInt k=0; k<LTraceIndexSets[1].rows(); k++ )
        rhsTemp[LTraceIndexSets[1][k]] = equivSrcs[1][k];

      // solve the system
      mCells[index]->solve( rhsTemp, resTemp );
      // add the sweep contribution to the result
      res[index] += resTemp;
      for( UInt k=0; k<indexSets1[5].size(); k++ )
        resTemp[indexSets1[5][k]] = tracesLR[3][index][k];
      for( UInt k=0; k<indexSets2[1].size(); k++ )
        resTemp[indexSets2[1][k]] = tracesUD[0][index][k];

     // use the solution to define the local trace DOFs that are required
     // we do need to store both the UpDown and LeftRight trace information
     // because both are needed for the L-sweep.
     // However, we only need to store the left and top traces because
     // those are the only ones needed in the next steps of this L sweep
      for( UInt k=0; k<indexSets1[1].rows(); k++ )
        tracesLR[0][index][k] += resTemp[indexSets1[1][k]];
      for( UInt k=0; k<indexSets1[2].rows(); k++ )
        tracesLR[1][index][k] += resTemp[indexSets1[2][k]];
      for( UInt k=0; k<indexSets2[4].rows(); k++ )
        tracesUD[2][index][k] += resTemp[indexSets2[4][k]];
      for( UInt k=0; k<indexSets2[5].rows(); k++ )
        tracesUD[3][index][k] += resTemp[indexSets2[5][k]];
    }

    // if we are not in the first column or first row and the element on the left
    // is not in the same process, we need to send the trace info to the left
    if( JJ>0 && mIndToProc(II,JJ)!=mIndToProc(II,JJ-1) && II>0 )
    {
      // devide up the trace0 info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR0_re(tracesLR[0][index].size());
      DoubleVector tracesLR0_im(tracesLR[0][index].size());
      for( UInt k=0; k<tracesLR0_re.size(); k++ )
      {
        tracesLR0_re[k] = real(tracesLR[0][index][k]);
        tracesLR0_im[k] = imag(tracesLR[0][index][k]);
      }
      // send the info
      MPI_Send( tracesLR0_re.data(), tracesLR0_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+0, mComm );
      MPI_Send( tracesLR0_im.data(), tracesLR0_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+1, mComm );

      // devide up the traceNP info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR1_re(tracesLR[1][index].size());
      DoubleVector tracesLR1_im(tracesLR[1][index].size());
      for( UInt k=0; k<tracesLR1_re.size(); k++ )
      {
        tracesLR1_re[k] = real(tracesLR[1][index][k]);
        tracesLR1_im[k] = imag(tracesLR[1][index][k]);
      }
      // send the info
      MPI_Send( tracesLR1_re.data(), tracesLR1_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+2, mComm );
      MPI_Send( tracesLR1_im.data(), tracesLR1_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+3, mComm );
    }

    // if we are not in the last row or last column and the element above
    // is not in the same process, we need to send the trace info to the top
    if( II<mNCells1-1 && mIndToProc(II,JJ)!=mIndToProc(II+1,JJ) && JJ<mNCells2-1 )
    {
      // devide up the traceN info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD2_re(tracesUD[2][index].size());
      DoubleVector tracesUD2_im(tracesUD[2][index].size());
      for( UInt k=0; k<tracesUD2_re.size(); k++ )
      {
        tracesUD2_re[k] = real(tracesUD[2][index][k]);
        tracesUD2_im[k] = imag(tracesUD[2][index][k]);
      }
      // send the info
      MPI_Send( tracesUD2_re.data(), tracesUD2_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+0, mComm );
      MPI_Send( tracesUD2_im.data(), tracesUD2_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+1, mComm );

      // devide up the traceNP info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD3_re(tracesUD[3][index].size());
      DoubleVector tracesUD3_im(tracesUD[3][index].size());
      for( UInt k=0; k<tracesUD3_re.size(); k++ )
      {
        tracesUD3_re[k] = real(tracesUD[3][index][k]);
        tracesUD3_im[k] = imag(tracesUD[3][index][k]);
      }
      // send the info
      MPI_Send( tracesUD3_re.data(), tracesUD3_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+2, mComm );
      MPI_Send( tracesUD3_im.data(), tracesUD3_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+3, mComm );
    }
  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}

// the functiont that adds the contribution of the diagonal sweep from bottom left to top right
void DDCheckerboard::addBL2TRSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> tracesLR, vector<vector<ComplexVector>> tracesUD, bool shiftedLR, bool shiftedUD ) const
{
  // get the local rank
  UInt rank;
  MPI_Comm_rank(mComm, (int*)&rank);
  
  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    diagInds[i] = I+J;
    locInds(I,J)=i;
  }
 
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 

  // get the number of cells in each direction
  UInt nCells1 = mNCells1; 
  UInt nCells2 = mNCells2; 

  // now loop over all local subdomains in the order they should be processed
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the local index of the subdomain to be processed and its row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
   
    // only do something if we are not in the first row or the first column
    // (these cases were handled by the right and up sweep)
    if( II>0 && JJ>0 )
    {
      // get the index sets for the local subdomain
      vector<UIntVector> indexSets1, indexSets2;
      mCells[index]->getIndexSets("LeftRight",shiftedLR,indexSets1);
      mCells[index]->getIndexSets("UpDown",shiftedUD,indexSets2);

      // get the trace information from the left and bottom element
      //
      // if the left element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ-1) )
      {
        UInt index_m1 = locInds(II,JJ-1);
        tracesLR[0][index] = tracesLR[2][index_m1];
        tracesLR[1][index] = tracesLR[3][index_m1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesLR0_re(indexSets1[1].size());   
        DoubleVector tracesLR0_im(indexSets1[1].size());   
        DoubleVector tracesLR1_re(indexSets1[2].size());   
        DoubleVector tracesLR1_im(indexSets1[2].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesLR0_re.data(), indexSets1[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesLR0_im.data(), indexSets1[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesLR1_re.data(), indexSets1[2].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesLR1_im.data(), indexSets1[2].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesLR[0][index] = tracesLR0_re+im*tracesLR0_im;
        tracesLR[1][index] = tracesLR1_re+im*tracesLR1_im;
      }
      // if the bottom element is in the same process, we just have to copy 
      // the information over
      if( mIndToProc(II,JJ)==mIndToProc(II-1,JJ) )
      {
        UInt index_m1 = locInds(II-1,JJ);
        tracesUD[0][index] = tracesUD[2][index_m1];
        tracesUD[1][index] = tracesUD[3][index_m1];
      }
      // if not, then we have to receive them from the correct process
      else
      {
        // alocate memory 
        DoubleVector tracesUD0_re(indexSets2[1].size());   
        DoubleVector tracesUD0_im(indexSets2[1].size());   
        DoubleVector tracesUD1_re(indexSets2[2].size());   
        DoubleVector tracesUD1_im(indexSets2[2].size());  
        // receive the traces and store them into local variables
        MPI_Status status;
        MPI_Recv( tracesUD0_re.data(), indexSets2[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesUD0_im.data(), indexSets2[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesUD1_re.data(), indexSets2[2].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesUD1_im.data(), indexSets2[2].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+3, mComm, &status );
        // use the local variables to define the traces
        Complex im(0.0,1.0);
        tracesUD[0][index] = tracesUD0_re+im*tracesUD0_im;
        tracesUD[1][index] = tracesUD1_re+im*tracesUD1_im;
      }
    
      // extract the L-shaped traces
      vector<ComplexVector> LTraceVals;
      mCells[index]->extractLTraces( "BottomLeft", shiftedLR, shiftedUD, tracesLR[0][index], tracesLR[1][index], tracesUD[0][index], tracesUD[1][index], LTraceVals );

      // compute the equivalent sources on the boundary (the right-hand side contributions to compute the polarized wave field)
      vector<ComplexVector> equivSrcs;
      mCells[index]->getEquivalentSources( "BottomLeft", shiftedLR, shiftedUD, LTraceVals[0], LTraceVals[1], equivSrcs );

      // get the indicex sets corresponding to L-shaped traces
      vector<UIntVector> LTraceIndexSets;
      mCells[index]->getLTraceIndexSets( "BottomLeft", shiftedLR, shiftedUD, LTraceIndexSets );

      // now we can solve the local system using the equivalent sources
      // set the rhs and allocate memory for the result
      ComplexVector rhsTemp(mCells[index]->getSize());
      for( UInt i=0; i<rhsTemp.size(); i++ )
        rhsTemp[i] = 0;
      ComplexVector resTemp = rhsTemp;
      for( UInt k=0; k<LTraceIndexSets[0].rows(); k++ )
        rhsTemp[LTraceIndexSets[0][k]] = equivSrcs[0][k];
      for( UInt k=0; k<LTraceIndexSets[1].rows(); k++ )
        rhsTemp[LTraceIndexSets[1][k]] = equivSrcs[1][k];

      // solve the system
      mCells[index]->solve( rhsTemp, resTemp );
      for( UInt k=0; k<indexSets1[1].size(); k++ )
        resTemp[indexSets1[1][k]] = tracesLR[0][index][k];
      for( UInt k=0; k<indexSets2[1].size(); k++ )
        resTemp[indexSets2[1][k]] = tracesUD[0][index][k];

      // add the sweep contribution to the result
      res[index] += resTemp;

     // use the solution to define the local trace DOFs that are required
     // we do need to store both the UpDown and LeftRight trace information
     // because both are needed for the L-sweep.
     // However, we only need to store the right and top traces because
     // those are the only ones needed in the next steps of this L sweep
      for( UInt k=0; k<indexSets1[4].rows(); k++ )
        tracesLR[2][index][k] += resTemp[indexSets1[4][k]];
      for( UInt k=0; k<indexSets1[5].rows(); k++ )
        tracesLR[3][index][k] += resTemp[indexSets1[5][k]];
      for( UInt k=0; k<indexSets2[4].rows(); k++ )
        tracesUD[2][index][k] += resTemp[indexSets2[4][k]];
      for( UInt k=0; k<indexSets2[5].rows(); k++ )
        tracesUD[3][index][k] += resTemp[indexSets2[5][k]];
    }

    // if we are not in the last column or first row and the element on the right
    // is not in the same process, we need to send the trace info to the right
    if( JJ<mNCells2-1 && mIndToProc(II,JJ)!=mIndToProc(II,JJ+1) && II>0 )
    {
      // devide up the traceN info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR2_re(tracesLR[2][index].size());
      DoubleVector tracesLR2_im(tracesLR[2][index].size());
      for( UInt k=0; k<tracesLR2_re.size(); k++ )
      {
        tracesLR2_re[k] = real(tracesLR[2][index][k]);
        tracesLR2_im[k] = imag(tracesLR[2][index][k]);
      }
      // send the info
      MPI_Send( tracesLR2_re.data(), tracesLR2_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+0, mComm );
      MPI_Send( tracesLR2_im.data(), tracesLR2_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+1, mComm );

      // devide up the traceNP info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesLR3_re(tracesLR[3][index].size());
      DoubleVector tracesLR3_im(tracesLR[3][index].size());
      for( UInt k=0; k<tracesLR3_re.size(); k++ )
      {
        tracesLR3_re[k] = real(tracesLR[3][index][k]);
        tracesLR3_im[k] = imag(tracesLR[3][index][k]);
      }
      // send the info
      MPI_Send( tracesLR3_re.data(), tracesLR3_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+2, mComm );
      MPI_Send( tracesLR3_im.data(), tracesLR3_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+3, mComm );
    }

    // if we are not in the last row or first column and the element above
    // is not in the same process, we need to send the trace info to the top
    if( II<mNCells1-1 && mIndToProc(II,JJ)!=mIndToProc(II+1,JJ) && JJ>0 )
    {
      // devide up the traceN info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD2_re(tracesUD[2][index].size());
      DoubleVector tracesUD2_im(tracesUD[2][index].size());
      for( UInt k=0; k<tracesUD2_re.size(); k++ )
      {
        tracesUD2_re[k] = real(tracesUD[2][index][k]);
        tracesUD2_im[k] = imag(tracesUD[2][index][k]);
      }
      // send the info
      MPI_Send( tracesUD2_re.data(), tracesUD2_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+0, mComm );
      MPI_Send( tracesUD2_im.data(), tracesUD2_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+1, mComm );

      // devide up the traceNP info into real and imaginary part to prepare
      // for sending
      DoubleVector tracesUD3_re(tracesUD[3][index].size());
      DoubleVector tracesUD3_im(tracesUD[3][index].size());
      for( UInt k=0; k<tracesUD3_re.size(); k++ )
      {
        tracesUD3_re[k] = real(tracesUD[3][index][k]);
        tracesUD3_im[k] = imag(tracesUD[3][index][k]);
      }
      // send the info
      MPI_Send( tracesUD3_re.data(), tracesUD3_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+2, mComm );
      MPI_Send( tracesUD3_im.data(), tracesUD3_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+3, mComm );
    }
  } // loop over all local elements is done

  // wait until everybody is done with their loop
  MPI_Barrier(mComm);
}

// the function to apply the system matrix
void DDCheckerboard::applySystemMatrix( const vector<ComplexVector> &x, vector<ComplexVector> &res ) const
{
  // first expand the input vector
  vector<ComplexVector> xEx;
  expand( x, xEx );

  // communicate the trace information to the neighboring elements
  // these DOFs are needed, since the stencil in the interior of the 
  // domain reaches into those DOFs. Once these are communicated the
  // system matrices can be applied locally and will result in the correct
  // values for the interior of the domain
  communicateBL2TRTraces( xEx, false );  
  communicateTR2BLTraces( xEx, false );  


  // apply the local system matrices
  vector<ComplexVector> resTemp(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    mCells[i]->applySystemMatrix( xEx[i], resTemp[i] ); 
  }

  // reduce the dofs back
  reduce(resTemp,res);
}

// the function that communicates trace information from bottom left to top right
void DDCheckerboard::communicateBL2TRTraces( vector<ComplexVector> &xEx, bool shifted ) const
{
  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    diagInds[i] = I+J;
    
    locInds(I,J)=i;
  }

  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 
  MPI_Barrier(mComm);

  // loop over the cells in the right order and exchange the information
  // of the vertical traces
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
   
    // get the trace indices
    vector<UIntVector> indexSets;
    mCells[index]->getIndexSets("LeftRight",shifted, indexSets );
    
    // if we are not in the first column, we have to receive the 
    // trace information from the left element
    if( JJ>0 )
    {
      // if the left element is in the same process, we just copy the information over
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ-1) )
      {
        // get the index of the left element
        UInt index_m1 = locInds(II,JJ-1);

        // get the trace indices of the left
        vector<UIntVector> indexSets_m1;
        mCells[index_m1]->getIndexSets("LeftRight",shifted, indexSets_m1 );
        // copy the information
        assert( indexSets[1].size()==indexSets_m1[4].size() );
        for( uint k=0; k<indexSets[1].size(); k++ )
        {
          xEx[index][indexSets[1][k]] = xEx[index_m1][indexSets_m1[4][k]];
        }
      }
      // if the left element is not in the same process, we have to receive info
      // from the process of the left element
      else
      {
        // allocate memory for the values to be received
        DoubleVector traces0_re( indexSets[1].size() );
        DoubleVector traces0_im( indexSets[1].size() );

        // receive the trace information from the left element
        MPI_Status status;
        MPI_Recv( traces0_re.data(), indexSets[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 2*JJ*mNCells1+2*II+0, mComm, &status );
        MPI_Recv( traces0_im.data(), indexSets[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 2*JJ*mNCells1+2*II+1, mComm, &status );

        // use the received information to set up the traces here
        Complex im(0.0,1.0);
        for( UInt k=0; k<indexSets[1].size(); k++ )
          xEx[index][indexSets[1][k]]=traces0_re[k]+im*traces0_im[k];
      }
    }

    // if we are not in the last column and the element on the right is not in the same process, we have to send the trace info on the
    // right side of the element
    if( JJ<mNCells2-1 )
    {
      if( mIndToProc(II,JJ)!=mIndToProc(II,JJ+1) )
      {
        // copy the trace info into buffer variables for sending
        DoubleVector traceN_re( indexSets[4].size() );
        DoubleVector traceN_im( indexSets[4].size() );
        for( UInt k=0; k<indexSets[4].size(); k++ )
        {
          traceN_re[k] = real(xEx[index][indexSets[4][k]]);
          traceN_im[k] = imag(xEx[index][indexSets[4][k]]);
        }

        // use the buffer variables to send
        MPI_Send( traceN_re.data(), traceN_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 2*(JJ+1)*mNCells1+2*II+0, mComm );
        MPI_Send( traceN_im.data(), traceN_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 2*(JJ+1)*mNCells1+2*II+1, mComm );
      }
    }
  }
  MPI_Barrier(mComm);

  // loop over the cells in the right order and exchange the information
  // of the horizontal traces
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
        
    // get the trace indices
    vector<UIntVector> indexSets;
    mCells[index]->getIndexSets("UpDown",shifted,indexSets );
   
    // if we are not in the first row, we have to receive the 
    // trace information from the bottom element
    if( II>0 )
    {
      // if the bottom element is in the same process, we just copy the information over
      if( mIndToProc(II,JJ)==mIndToProc(II-1,JJ) )
      {
        // get the index of the bottom element
        UInt index_m1 = locInds(II-1,JJ);

        // get the trace indices of the bottom element
        vector<UIntVector> indexSets_m1;
        mCells[index_m1]->getIndexSets("UpDown",shifted,indexSets_m1 );

        // copy the information
        assert( indexSets[1].size()==indexSets_m1[4].size() );
        for( uint k=0; k<indexSets[1].size(); k++ )
          xEx[index][indexSets[1][k]] = xEx[index_m1][indexSets_m1[4][k]];
      }
      // if the botto element is not in the same process, we have to receive it from the
      // process of the bottom element
      else
      {
        // allocate memory for the values to be received
        DoubleVector traces0_re( indexSets[1].size() );
        DoubleVector traces0_im( indexSets[1].size() );
 
        // receive the trace information from the bottom element
        MPI_Status status;
        MPI_Recv( traces0_re.data(), indexSets[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 2*JJ*mNCells1+2*II+0, mComm, &status );
        MPI_Recv( traces0_im.data(), indexSets[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 2*JJ*mNCells1+2*II+1, mComm, &status );

        // use the received information to set up the traces here
        Complex im(0.0,1.0);
        for( UInt k=0; k<indexSets[1].size(); k++ )
          xEx[index][indexSets[1][k]]=traces0_re[k]+im*traces0_im[k];
      }
    }

    // if we are not in the last row and the element on the top is not in the same process, we have to send the trace info on the
    // top side of the element
    if( II<mNCells1-1 )
    {
      if( mIndToProc(II,JJ)!=mIndToProc(II+1,JJ) )
      {
        // copy the trace info into buffer variables for sending
        DoubleVector traceN_re( indexSets[4].size() );
        DoubleVector traceN_im( indexSets[4].size() );
        for( UInt k=0; k<indexSets[4].size(); k++ )
        {
          traceN_re[k] = real(xEx[index][indexSets[4][k]]);
          traceN_im[k] = imag(xEx[index][indexSets[4][k]]);
        }
        
        // use the buffer variables to send
        MPI_Send( traceN_re.data(), traceN_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 2*JJ*mNCells1+2*(II+1)+0, mComm );
        MPI_Send( traceN_im.data(), traceN_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 2*JJ*mNCells1+2*(II+1)+1, mComm );
      }
    }
  }
  MPI_Barrier(mComm);
}

// the function that communicates trace information from top right to bottom left
void DDCheckerboard::communicateTR2BLTraces( vector<ComplexVector> &xEx, bool shifted ) const
{
  // first find out in what row and column the local cells are located
  // and compute their diagonal index
  // note that we define the diagonals from top right to bottom left
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    I = mNCells1-1-I;
    J = mNCells2-1-J;

    diagInds[i] = I+J;
    
    I = mNCells1-1-I;
    J = mNCells2-1-J;
    
    locInds(I,J)=i;
  }

  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 
  MPI_Barrier(mComm);

  // loop over the cells in the right order and exchange the information
  // of the vertical traces
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
        
    // get the trace indices
    vector<UIntVector> indexSets;
    mCells[index]->getIndexSets("LeftRight",shifted, indexSets );
   
    // if we are not in the last column, we have to receive the 
    // trace information from the right element
    if( JJ<mNCells2-1 )
    {
      // if the right element is in the same process, we just copy the info over
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ+1) )
      {
        // get the index of the right element
        UInt index_p1 = locInds(II,JJ+1);

        // get the trace indices of the right element
        vector<UIntVector> indexSets_p1;
        mCells[index_p1]->getIndexSets("LeftRight",shifted, indexSets_p1 );

        // copy the informaiton
        assert( indexSets[5].size()==indexSets_p1[2].size() );
        for( uint k=0; k<indexSets[5].size(); k++ )
          xEx[index][indexSets[5][k]] = xEx[index_p1][indexSets_p1[2][k]];
      }
      // if the right element is not in the same process, we have to receive info
      // from the process of the right element
      else
      {

        // allocate memory for the values to be received
        DoubleVector tracesNP_re( indexSets[5].size() );
        DoubleVector tracesNP_im( indexSets[5].size() );

        // receive the trace information from the right element
        MPI_Status status;
        MPI_Recv( tracesNP_re.data(), indexSets[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 2*JJ*mNCells1+2*II+0, mComm, &status );
        MPI_Recv( tracesNP_im.data(), indexSets[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 2*JJ*mNCells1+2*II+1, mComm, &status );

        // use the received information to set up the traces here
        Complex im(0.0,1.0);
        for( UInt k=0; k<indexSets[5].size(); k++ )
          xEx[index][indexSets[5][k]]=tracesNP_re[k]+im*tracesNP_im[k];
      }
    }

    // if we are not in the first column and the element on the left is not in the same
    // process, we have to send the trace info on the right side of the element
    if( JJ>0 )
    {
      if( mIndToProc(II,JJ)!=mIndToProc(II,JJ-1) )
      {
        // copy the trace info into buffer variables for sending
        DoubleVector trace1_re( indexSets[2].size() );
        DoubleVector trace1_im( indexSets[2].size() );
        for( UInt k=0; k<indexSets[2].size(); k++ )
        {
          trace1_re[k] = real(xEx[index][indexSets[2][k]]);
          trace1_im[k] = imag(xEx[index][indexSets[2][k]]);
        }

        // use the buffer variables to send
        MPI_Send( trace1_re.data(), trace1_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 2*(JJ-1)*mNCells1+2*II+0, mComm );
        MPI_Send( trace1_im.data(), trace1_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 2*(JJ-1)*mNCells1+2*II+1, mComm );
      }
    }
  }
  MPI_Barrier(mComm);

  // loop over the cells inthe right order and exchange th information of the
  // horizontal traces
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the row and column index
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
        
    // get the trace indices
    vector<UIntVector> indexSets;
    mCells[index]->getIndexSets("UpDown",shifted,indexSets );
  
    // if we are not in the last row, we have to receive
    // the trace information from the bottom element
    if( II<mNCells1-1 )
    {
      // if the top element is in the same process, we just copy the informaiton over
      if( mIndToProc(II,JJ)==mIndToProc(II+1,JJ) )
      {
        // get the index of the top element
        UInt index_p1 = locInds(II+1,JJ);

        // get the trace indices of the top element
        vector<UIntVector> indexSets_p1;
        mCells[index_p1]->getIndexSets("UpDown",shifted,indexSets_p1 );

        // copy the information
        assert( indexSets[5].size()==indexSets_p1[2].size() );
        for( uint k=0; k<indexSets[5].size(); k++ )
          xEx[index][indexSets[5][k]] = xEx[index_p1][indexSets_p1[2][k]];
      }
      // if the top element is not in the same process, we have to receive info from the
      // process of the top element
      else
      {
        // allocate memory for the trace values to be received
        DoubleVector tracesNP_re( indexSets[5].size() );
        DoubleVector tracesNP_im( indexSets[5].size() );

        // receive the trace information from the top element
        MPI_Status status;
        MPI_Recv( tracesNP_re.data(), indexSets[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 2*JJ*mNCells1+2*II+0, mComm, &status );
        MPI_Recv( tracesNP_im.data(), indexSets[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 2*JJ*mNCells1+2*II+1, mComm, &status );

        // use the received information to set up the traces here
        Complex im(0.0,1.0);
        for( UInt k=0; k<indexSets[5].size(); k++ )
          xEx[index][indexSets[5][k]]=tracesNP_re[k]+im*tracesNP_im[k];
      }
    }

    // if we are not in the first row and the element on the bottom is not in the same process, we have to send the trace info on the
    // bottom side of the element
    if( II>0 )
    {
      if( mIndToProc(II,JJ)!=mIndToProc(II-1,JJ) )
      {
        // copy the trace info into buffer variables for sending
        DoubleVector trace1_re( indexSets[2].size() );
        DoubleVector trace1_im( indexSets[2].size() );
        for( UInt k=0; k<indexSets[2].size(); k++ )
        {
          trace1_re[k] = real(xEx[index][indexSets[2][k]]);
          trace1_im[k] = imag(xEx[index][indexSets[2][k]]);
        }

        // use the buffer variables to send
        MPI_Send( trace1_re.data(), trace1_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 2*JJ*mNCells1+2*(II-1)+0, mComm );
        MPI_Send( trace1_im.data(), trace1_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 2*JJ*mNCells1+2*(II-1)+1, mComm );
      }
    }
  }
  MPI_Barrier(mComm);
}

// the function that converts non-shifted-DOF-storage to shifted-DOF-storage
void DDCheckerboard::convertNonShiftedToShiftedLR( vector<ComplexVector> &x, bool shiftedUD ) const
{
  // first expand the DOFs 
  vector<ComplexVector> xEx;
  expand( x, xEx, false, shiftedUD );
  
  // we need to transfer information from top right to bottom left
  // we now find the order of the elements to be processed
  // start by assiging diagonal indices to each element
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    I = mNCells1-1-I;
    J = mNCells2-1-J;

    diagInds[i] = I+J;
    
    I = mNCells1-1-I;
    J = mNCells2-1-J;
    
    locInds(I,J)=i;
  }
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 
  MPI_Barrier(mComm);

  // now loop over the subdomains in the right order and exchange the information
  // of the vertical traces
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the row and column indices
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
       
    // get the trace indices
    vector<UIntVector> indexSets;
    vector<UIntVector> indexSetsShifted;
    mCells[index]->getIndexSets("LeftRight",false, indexSets );
    mCells[index]->getIndexSets("LeftRight",true, indexSetsShifted );
   
    // if we are no the last column, we have to receive info from the right element
    if( JJ<mNCells2-1 )
    {
      // if the right element is in the same process, we just have to copy the info
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ+1) )
      {
        // get the index of the right element
        UInt index_p1 = locInds(II,JJ+1);

        // get the trace indices of the elment on the right
        vector<UIntVector> indexSets_p1;
        vector<UIntVector> indexSetsShifted_p1;
        mCells[index_p1]->getIndexSets("LeftRight",false, indexSets_p1 );
        mCells[index_p1]->getIndexSets("LeftRight",true, indexSetsShifted_p1 );

        // copy the info
        assert( indexSets[5].size()==indexSets_p1[2].size() );
        for( uint k=0; k<indexSets[5].size(); k++ )
          xEx[index][indexSets[5][k]] = xEx[index_p1][indexSets_p1[2][k]];
        assert( indexSetsShifted[4].size()==indexSetsShifted_p1[1].size() );
        for( uint k=0; k<indexSetsShifted[4].size(); k++ )
          xEx[index][indexSetsShifted[4][k]] = xEx[index_p1][indexSetsShifted_p1[1][k]];
      }
      // if the right element is not in the same process, we have to receive the info
      // from the element on the right
      else
      {
        // allocate memory for the trace info to be received
        DoubleVector tracesNP_re( indexSets[5].size() );
        DoubleVector tracesNP_im( indexSets[5].size() );
        DoubleVector tracesNShifted_re( indexSetsShifted[4].size() );
        DoubleVector tracesNShifted_im( indexSetsShifted[4].size() );
 
        // receive the info
        MPI_Status status;
        MPI_Recv( tracesNP_re.data(), indexSets[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesNP_im.data(), indexSets[5].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesNShifted_re.data(), indexSetsShifted[4].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesNShifted_im.data(), indexSetsShifted[4].size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*JJ*mNCells1+4*II+3, mComm, &status );

        // use the received information to set up the DOFs here
        Complex im(0.0,1.0);
        for( UInt k=0; k<indexSets[5].size(); k++ )
          xEx[index][indexSets[5][k]]=tracesNP_re[k]+im*tracesNP_im[k];
        for( UInt k=0; k<indexSetsShifted[4].size(); k++ )
          xEx[index][indexSetsShifted[4][k]]=tracesNShifted_re[k]+im*tracesNShifted_im[k];
      }
    }

    // if we are not in the first column and the element on the left is not
    // in the same process, we have to send the info there
    if( JJ>0 )
    {
      if( mIndToProc(II,JJ)!=mIndToProc(II,JJ-1) )
      {
        // copy the trace info to be sent into buffer variables
        DoubleVector trace1_re( indexSets[2].size() );
        DoubleVector trace1_im( indexSets[2].size() );
        DoubleVector trace0Shifted_re( indexSetsShifted[1].size() );
        DoubleVector trace0Shifted_im( indexSetsShifted[1].size() );
        for( UInt k=0; k<indexSets[2].size(); k++ )
        {
          trace1_re[k] = real(xEx[index][indexSets[2][k]]);
          trace1_im[k] = imag(xEx[index][indexSets[2][k]]);
        }
        for( UInt k=0; k<indexSetsShifted[1].size(); k++ )
        {
          trace0Shifted_re[k] = real(xEx[index][indexSetsShifted[1][k]]);
          trace0Shifted_im[k] = imag(xEx[index][indexSetsShifted[1][k]]);
        }

        // send the info
        MPI_Send( trace1_re.data(), trace1_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+0, mComm );
        MPI_Send( trace1_im.data(), trace1_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+1, mComm );
        MPI_Send( trace0Shifted_re.data(), trace0Shifted_re.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+2, mComm );
        MPI_Send( trace0Shifted_im.data(), trace0Shifted_im.size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*(JJ-1)*mNCells1+4*II+3, mComm );
      }
    }
  }
  MPI_Barrier(mComm);

  // reduce the dofs to the shifted ones
  reduce(xEx,x,true,shiftedUD);
}

// the function that converts non-shifted-DOF-storage to shifted-DOF-storage
void DDCheckerboard::convertNonShiftedToShiftedUD( vector<ComplexVector> &x, bool shiftedLR ) const
{
  // first expand the DOFs 
  vector<ComplexVector> xEx;
  expand( x, xEx, shiftedLR, false );
  
  // we need to transfer information from top right to bottom left
  // we now find the order of the elements to be processed
  // start by assiging diagonal indices to each element
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    I = mNCells1-1-I;
    J = mNCells2-1-J;

    diagInds[i] = I+J;
    
    I = mNCells1-1-I;
    J = mNCells2-1-J;
    
    locInds(I,J)=i;
  }
  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 
  MPI_Barrier(mComm);

  // now loop over the subdomains in the right order and exchange the information
  // of the horizontal traces
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the row and column indices
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
        
    // get the trace indices
    vector<UIntVector> indexSets;
    vector<UIntVector> indexSetsShifted;
    mCells[index]->getIndexSets("UpDown",false,indexSets );
    mCells[index]->getIndexSets("UpDown",true,indexSetsShifted );
   
    // if we are no the last row, we have to receive info from the top element
    if( II<mNCells1-1 )
    {
      // if the top element is in the same process, we just have to copy the info
      if( mIndToProc(II,JJ)==mIndToProc(II+1,JJ) )
      {
        // get the index of the top element
        UInt index_p1 = locInds(II+1,JJ);

        // get the trace indices of the top element
        vector<UIntVector> indexSets_p1;
        vector<UIntVector> indexSetsShifted_p1;
        mCells[index_p1]->getIndexSets("UpDown",false,indexSets_p1 );
        mCells[index_p1]->getIndexSets("UpDown",true,indexSetsShifted_p1 );

        // copy the info
        assert( indexSets[5].size()==indexSets_p1[2].size() );
        for( uint k=0; k<indexSets[5].size(); k++ )
          xEx[index][indexSets[5][k]] = xEx[index_p1][indexSets_p1[2][k]];
        assert( indexSetsShifted[4].size()==indexSetsShifted_p1[1].size() );
        for( uint k=0; k<indexSetsShifted[4].size(); k++ )
          xEx[index][indexSetsShifted[4][k]] = xEx[index_p1][indexSetsShifted_p1[1][k]];
      }
      // if the top element is not in the same process, we have to receive the info
      // from the element on the top
      else
      {
        // allocate memory for the trace info to be received
        DoubleVector tracesNP_re( indexSets[5].size() );
        DoubleVector tracesNP_im( indexSets[5].size() );
        DoubleVector tracesNShifted_re( indexSetsShifted[4].size() );
        DoubleVector tracesNShifted_im( indexSetsShifted[4].size() );
 
        // receive info
        MPI_Status status;
        MPI_Recv( tracesNP_re.data(), indexSets[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( tracesNP_im.data(), indexSets[5].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( tracesNShifted_re.data(), indexSetsShifted[4].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( tracesNShifted_im.data(), indexSetsShifted[4].size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*II+3, mComm, &status );

        // use the received information to set up the DOFs here
        Complex im(0.0,1.0);
        for( UInt k=0; k<indexSets[5].size(); k++ )
          xEx[index][indexSets[5][k]]=tracesNP_re[k]+im*tracesNP_im[k];
        for( UInt k=0; k<indexSetsShifted[4].size(); k++ )
          xEx[index][indexSetsShifted[4][k]]=tracesNShifted_re[k]+im*tracesNShifted_im[k];
      }
    }

    // if we are not in the first row and the element on the top is not
    // in the same process, we have to send the info there
    if( II>0 )
    {
      if( mIndToProc(II,JJ)!=mIndToProc(II-1,JJ) )
      {
        // copy the trace info to be sent into buffer variables
        DoubleVector trace1_re( indexSets[2].size() );
        DoubleVector trace1_im( indexSets[2].size() );
        DoubleVector trace0Shifted_re( indexSetsShifted[1].size() );
        DoubleVector trace0Shifted_im( indexSetsShifted[1].size() );
        for( UInt k=0; k<indexSets[2].size(); k++ )
        {
          trace1_re[k] = real(xEx[index][indexSets[2][k]]);
          trace1_im[k] = imag(xEx[index][indexSets[2][k]]);
        }
        for( UInt k=0; k<indexSetsShifted[1].size(); k++ )
        {
          trace0Shifted_re[k] = real(xEx[index][indexSetsShifted[1][k]]);
          trace0Shifted_im[k] = imag(xEx[index][indexSetsShifted[1][k]]);
        }

        // send info
        MPI_Send( trace1_re.data(), trace1_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+0, mComm );
        MPI_Send( trace1_im.data(), trace1_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+1, mComm );
        MPI_Send( trace0Shifted_re.data(), trace0Shifted_re.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+2, mComm );
        MPI_Send( trace0Shifted_im.data(), trace0Shifted_im.size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*(II-1)+3, mComm );
      }
    }
  }
  MPI_Barrier(mComm);

  // reduce the dofs to the shifted ones
  reduce(xEx,x,shiftedLR,true);
}

// the function that converts shifted-DOF-storage to non-shifted-DOF-storage
void DDCheckerboard::convertShiftedToNonShiftedLR( vector<ComplexVector> &x, bool shiftedUD ) const
{
  // first expand the DOFs 
  vector<ComplexVector> xEx;
  expand( x, xEx, true, shiftedUD );
  
  // we need to transfer information from bottom left to top right
  // we now find the order of the elements to be processed
  // start by assiging diagonal indices to each element
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    diagInds[i] = I+J;
    
    locInds(I,J)=i;
  }

  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 
  MPI_Barrier(mComm);

  // now loop over the subdomains in the right order and exchange the information
  // of the vertical traces
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the row and column indices
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
        
    // get the trace indices
    vector<UIntVector> indexSets;
    vector<UIntVector> indexSetsShifted;
    mCells[index]->getIndexSets("LeftRight",false, indexSets );
    mCells[index]->getIndexSets("LeftRight",true, indexSetsShifted );
   
    // if we are not in the first column, we have to receive info from the left element
    if( JJ>0 )
    {
      // if the left element is in the same process, we just have to copy the info
      if( mIndToProc(II,JJ)==mIndToProc(II,JJ-1) )
      {
        // get the index of the left element
        UInt index_m1 = locInds(II,JJ-1);

        // get the trace indices of the left element
        vector<UIntVector> indexSets_m1;
        vector<UIntVector> indexSetsShifted_m1;
        mCells[index_m1]->getIndexSets("LeftRight",false, indexSets_m1 );
        mCells[index_m1]->getIndexSets("LeftRight",true, indexSetsShifted_m1 );

        // copy the info
        assert( indexSetsShifted[1].size()==indexSetsShifted_m1[4].size() );
        for( uint k=0; k<indexSetsShifted[1].size(); k++ )
          xEx[index][indexSetsShifted[1][k]] = xEx[index_m1][indexSetsShifted_m1[4][k]];
        assert( indexSets[2].size()==indexSets_m1[5].size() );
        for( uint k=0; k<indexSets[2].size(); k++ )
          xEx[index][indexSets[2][k]] = xEx[index_m1][indexSets_m1[5][k]];
      }
      // if the left element is not in the same process, we have to receive the info
      // from the element on the left
      else
      {
        // allocate memory for the trace info to be received
        DoubleVector traces0Shifted_re( indexSetsShifted[1].size() );
        DoubleVector traces0Shifted_im( indexSetsShifted[1].size() );
        DoubleVector traces1_re( indexSets[2].size() );
        DoubleVector traces1_im( indexSets[2].size() );
 
        // receive the info
        MPI_Status status;
        MPI_Recv( traces0Shifted_re.data(), indexSetsShifted[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( traces0Shifted_im.data(), indexSetsShifted[1].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( traces1_re.data(), indexSets[2].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( traces1_im.data(), indexSets[2].size(), MPI_DOUBLE, mIndToProc(II,JJ-1), 4*JJ*mNCells1+4*II+3, mComm, &status );

        // use the received information to set up the DOFs here
        Complex im(0.0,1.0);
        for( UInt k=0; k<indexSetsShifted[1].size(); k++ )
          xEx[index][indexSetsShifted[1][k]]=traces0Shifted_re[k]+im*traces0Shifted_im[k];
        for( UInt k=0; k<indexSets[2].size(); k++ )
          xEx[index][indexSets[2][k]]=traces1_re[k]+im*traces1_im[k];
      }
    }

    // if we are not in the last column and the element on the right is not
    // in the same process, we have to send the info there
    if( JJ<mNCells2-1 )
    {
      if( mIndToProc(II,JJ)!=mIndToProc(II,JJ+1) )
      {
        // copy the trace info to be sent into buffer variables
        DoubleVector traceNShifted_re( indexSetsShifted[4].size() );
        DoubleVector traceNShifted_im( indexSetsShifted[4].size() );
        DoubleVector traceNP_re( indexSets[5].size() );
        DoubleVector traceNP_im( indexSets[5].size() );
        for( UInt k=0; k<indexSetsShifted[4].size(); k++ )
        {
          traceNShifted_re[k] = real(xEx[index][indexSetsShifted[4][k]]);
          traceNShifted_im[k] = imag(xEx[index][indexSetsShifted[4][k]]);
        }
        for( UInt k=0; k<indexSets[5].size(); k++ )
        {
          traceNP_re[k] = real(xEx[index][indexSets[5][k]]);
          traceNP_im[k] = imag(xEx[index][indexSets[5][k]]);
        }

        // send the info
        MPI_Send( traceNShifted_re.data(), traceNShifted_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+0, mComm );
        MPI_Send( traceNShifted_im.data(), traceNShifted_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+1, mComm );
        MPI_Send( traceNP_re.data(), traceNP_re.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+2, mComm );
        MPI_Send( traceNP_im.data(), traceNP_im.size(), MPI_DOUBLE, mIndToProc(II,JJ+1), 4*(JJ+1)*mNCells1+4*II+3, mComm );
      }
    }
  }
  MPI_Barrier(mComm);

  // reduce the dofs to the non-shifted ones
  reduce(xEx,x,false,shiftedUD);
}

// the function that converts shifted-DOF-storage to non-shifted-DOF-storage
void DDCheckerboard::convertShiftedToNonShiftedUD( vector<ComplexVector> &x, bool shiftedLR ) const
{
  // first expand the DOFs 
  vector<ComplexVector> xEx;
  expand( x, xEx, shiftedLR, true );
  
  // we need to transfer information from bottom left to top right
  // we now find the order of the elements to be processed
  // start by assiging diagonal indices to each element
  UIntMatrix locInds(mNCells1,mNCells2);
  DoubleVector diagInds(mCells.size());
  for( UInt i=0; i<mCells.size(); i++ )
  {
    UInt I = mCells[i]->getRowIndex(); 
    UInt J = mCells[i]->getColIndex(); 

    diagInds[i] = I+J;
    
    locInds(I,J)=i;
  }

  // now use the diagonal index to find out in what order the local
  // cells should be processed
  UIntVector subdomainOrder;
  getSortedIndices( diagInds, subdomainOrder ); 
  MPI_Barrier(mComm);

  // now loop over the subdomains in the right order and exchange the information
  // of the horizontal traces
  for( UInt d=0; d<subdomainOrder.size(); d++ )
  {
    // get the row and column indices
    UInt index = subdomainOrder[d];
    UInt II = mCells[index]->getRowIndex();
    UInt JJ = mCells[index]->getColIndex();
        
    // get the trace indices
    vector<UIntVector> indexSets;
    vector<UIntVector> indexSetsShifted;
    mCells[index]->getIndexSets("UpDown",false,indexSets );
    mCells[index]->getIndexSets("UpDown",true,indexSetsShifted );
   
    // if we are not in the first row, we have to receive info from the bottom element
    if( II>0 )
    {
      // if the bottom element is in the same process, we just have to copy the info
      if( mIndToProc(II,JJ)==mIndToProc(II-1,JJ) )
      {
        // get the index of the bottom element
        UInt index_m1 = locInds(II-1,JJ);

        // get the trace indices of the bottom element
        vector<UIntVector> indexSets_m1;
        vector<UIntVector> indexSetsShifted_m1;
        mCells[index_m1]->getIndexSets("UpDown",false,indexSets_m1 );
        mCells[index_m1]->getIndexSets("UpDown",true,indexSetsShifted_m1 );

        // copy the info
        assert( indexSetsShifted[1].size()==indexSetsShifted_m1[4].size() );
        for( uint k=0; k<indexSetsShifted[1].size(); k++ )
          xEx[index][indexSetsShifted[1][k]] = xEx[index_m1][indexSetsShifted_m1[4][k]];
        assert( indexSets[2].size()==indexSets_m1[5].size() );
        for( uint k=0; k<indexSets[2].size(); k++ )
          xEx[index][indexSets[2][k]] = xEx[index_m1][indexSets_m1[5][k]];
      }
      // if the bottom element is not in the same process, we have to receive the info
      // from the element on the bottom
      else
      {
        // allocate memory for the trace info to be received
        DoubleVector traces0Shifted_re( indexSetsShifted[1].size() );
        DoubleVector traces0Shifted_im( indexSetsShifted[1].size() );
        DoubleVector traces1_re( indexSets[2].size() );
        DoubleVector traces1_im( indexSets[2].size() );
 
        // receive info
        MPI_Status status;
        MPI_Recv( traces0Shifted_re.data(), indexSetsShifted[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+0, mComm, &status );
        MPI_Recv( traces0Shifted_im.data(), indexSetsShifted[1].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+1, mComm, &status );
        MPI_Recv( traces1_re.data(), indexSets[2].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+2, mComm, &status );
        MPI_Recv( traces1_im.data(), indexSets[2].size(), MPI_DOUBLE, mIndToProc(II-1,JJ), 4*JJ*mNCells1+4*II+3, mComm, &status );

        // use the received information to set up the DOFs here
        Complex im(0.0,1.0);
        for( UInt k=0; k<indexSetsShifted[1].size(); k++ )
          xEx[index][indexSetsShifted[1][k]]=traces0Shifted_re[k]+im*traces0Shifted_im[k];
        for( UInt k=0; k<indexSets[2].size(); k++ )
          xEx[index][indexSets[2][k]]=traces1_re[k]+im*traces1_im[k];
      }
    }

    // if we are not in the last row and the element on the bottom is not
    // in the same process, we have to send the info there
    if( II<mNCells1-1 )
    {
      if( mIndToProc(II,JJ)!=mIndToProc(II+1,JJ) )
      {
        // copy the trace info to be sent into buffer variables
        DoubleVector traceNShifted_re( indexSetsShifted[4].size() );
        DoubleVector traceNShifted_im( indexSetsShifted[4].size() );
        DoubleVector traceNP_re( indexSets[5].size() );
        DoubleVector traceNP_im( indexSets[5].size() );
        for( UInt k=0; k<indexSetsShifted[4].size(); k++ )
        {
          traceNShifted_re[k] = real(xEx[index][indexSetsShifted[4][k]]);
          traceNShifted_im[k] = imag(xEx[index][indexSetsShifted[4][k]]);
        }
        for( UInt k=0; k<indexSets[5].size(); k++ )
        {
          traceNP_re[k] = real(xEx[index][indexSets[5][k]]);
          traceNP_im[k] = imag(xEx[index][indexSets[5][k]]);
        }
        
        // send info
        MPI_Send( traceNShifted_re.data(), traceNShifted_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+0, mComm );
        MPI_Send( traceNShifted_im.data(), traceNShifted_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+1, mComm );
        MPI_Send( traceNP_re.data(), traceNP_re.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+2, mComm );
        MPI_Send( traceNP_im.data(), traceNP_im.size(), MPI_DOUBLE, mIndToProc(II+1,JJ), 4*JJ*mNCells1+4*(II+1)+3, mComm );
      }
    }

  }
  MPI_Barrier(mComm);

  // reduce the dofs to the non-shifted ones
  reduce(xEx,x,shiftedLR,false);
}

// the function that returns the dimension of the underlying problem
UInt DDCheckerboard::getDim() const
{
  return mCells[0]->getDim();
}

// the function that returns the number of cells
UInt DDCheckerboard::getNCells() const
{
  return mCells.size();
}

// the function that returns a pointer to a subdomain (cell)
Subdomain* DDCheckerboard::getCell( UInt i ) const
{
  return mCells[i];
}

// the function that returns the global mesh size
Double DDCheckerboard::getMeshSize() const
{
  // get the curretn rank and total number of processors
  UInt nprocs,rank;
  MPI_Comm_size(mComm, (int*)&nprocs);
  MPI_Comm_rank(mComm, (int*)&rank);

  // get the minimum mesh size of all the cells stored in this rank
  Double h = mCells[0]->getMeshSize();
  for( UInt i=0; i<mCells.size(); i++ )
    h = min( h, mCells[i]->getMeshSize() );

  // gather all local mesh sizes in all ranks
  vector<Double> hs(nprocs);
  MPI_Allgather( &h, 1, MPI_DOUBLE, hs.data(), 1, MPI_DOUBLE, mComm );

  // determine the minimum of all mesh sizes
  for( UInt i=0; i<nprocs; i++ )
    h = min( h, hs[i] );

  // wait until everybody is done
  MPI_Barrier(mComm);

  // return the mesh size
  return h;
}

// compute the windowed right-hand side vector
void DDCheckerboard::getSeparatedRHS(const vector<ComplexVector> &rhs, vector<ComplexVector> &res, bool shiftedLR, bool shiftedUD ) const
{
  // expand the rhs and store it in the result
  expand(rhs,res);

  // depending on the shift of the traces set the corresponding dofs to zero
  // both traces are shifted
  if( !shiftedLR && !shiftedUD )
  {
    // loop over the subdomains
    for( UInt i=0; i<res.size(); i++ )
    {
      // set both non-shifted Left-Right traces to zero
      vector<UIntVector> indexSets;
      mCells[i]->getIndexSets("LeftRight",false,indexSets);
      for( UInt j=1; j<=2; j++ )
        for( UInt k=0; k<indexSets[j].size(); k++ )
          res[i][indexSets[j][k]]=0.0;
      for( UInt j=4; j<=5; j++ )
        for( UInt k=0; k<indexSets[j].size(); k++ )
          res[i][indexSets[j][k]]=0.0;
      
      // set both non-shifted Up-Down traces to zero
      mCells[i]->getIndexSets("UpDown",false,indexSets);
      for( UInt j=1; j<=2; j++ )
        for( UInt k=0; k<indexSets[j].size(); k++ )
          res[i][indexSets[j][k]]=0.0;
      for( UInt j=4; j<=5; j++ )
        for( UInt k=0; k<indexSets[j].size(); k++ )
          res[i][indexSets[j][k]]=0.0;
    }
  }
  // the Left-Right trace is shifted
  else if( !shiftedLR )
  {
    // loop over the subdomains
    for( UInt i=0; i<res.size(); i++ )
    {
      // set both non-shifted Left-Right traces to zero
      vector<UIntVector> indexSets;
      mCells[i]->getIndexSets("LeftRight",false,indexSets);
      for( UInt j=1; j<=2; j++ )
        for( UInt k=0; k<indexSets[j].size(); k++ )
          res[i][indexSets[j][k]]=0.0;
      for( UInt j=4; j<=5; j++ )
        for( UInt k=0; k<indexSets[j].size(); k++ )
          res[i][indexSets[j][k]]=0.0;
      
      // set all dofs to zero that are not the non-shifted Up-Down traces
      mCells[i]->getIndexSets("UpDown",false,indexSets);
      for( UInt k=0; k<indexSets[0].size(); k++ )
        res[i][indexSets[0][k]]=0.0;
      for( UInt k=0; k<indexSets[3].size(); k++ )
        res[i][indexSets[3][k]]=0.0;
      for( UInt k=0; k<indexSets[6].size(); k++ )
        res[i][indexSets[6][k]]=0.0;
    }
  }
  // the Up-Down trace is shifted
  else if( !shiftedUD )
  {
    // loop over the subdomains
    for( UInt i=0; i<res.size(); i++ )
    {
      // set all dofs to zero that are not the non-shifted Left-Right traces
      vector<UIntVector> indexSets;
      mCells[i]->getIndexSets("LeftRight",false,indexSets);
      for( UInt k=0; k<indexSets[0].size(); k++ )
        res[i][indexSets[0][k]]=0.0;
      for( UInt k=0; k<indexSets[3].size(); k++ )
        res[i][indexSets[3][k]]=0.0;
      for( UInt k=0; k<indexSets[6].size(); k++ )
        res[i][indexSets[6][k]]=0.0;
      
      // set both non-shifted Up-Down traces to zero
      mCells[i]->getIndexSets("UpDown",false,indexSets);
      for( UInt j=1; j<=2; j++ )
        for( UInt k=0; k<indexSets[j].size(); k++ )
          res[i][indexSets[j][k]]=0.0;
      for( UInt j=4; j<=5; j++ )
        for( UInt k=0; k<indexSets[j].size(); k++ )
          res[i][indexSets[j][k]]=0.0;
    }
  }
  // no traces are shifted
  else
  {
    // loop over the subdomains
    for( UInt i=0; i<res.size(); i++ )
    {
      // set all dofs to zero that are not the non-shifted Left-Right traces
      vector<UIntVector> indexSets;
      mCells[i]->getIndexSets("LeftRight",false,indexSets);
      for( UInt k=0; k<indexSets[0].size(); k++ )
        res[i][indexSets[0][k]]=0.0;
      for( UInt k=0; k<indexSets[3].size(); k++ )
        res[i][indexSets[3][k]]=0.0;
      for( UInt k=0; k<indexSets[6].size(); k++ )
        res[i][indexSets[6][k]]=0.0;
      
      // set all dofs to zero that are not the non-shifted Up-Down traces
      mCells[i]->getIndexSets("UpDown",false,indexSets);
      for( UInt k=0; k<indexSets[0].size(); k++ )
        res[i][indexSets[0][k]]=0.0;
      for( UInt k=0; k<indexSets[3].size(); k++ )
        res[i][indexSets[3][k]]=0.0;
      for( UInt k=0; k<indexSets[6].size(); k++ )
        res[i][indexSets[6][k]]=0.0;
    }
  }
}

