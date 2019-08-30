#include "Input.h"

#include <iostream>
#include "libconfig.h++"
#include "utils.h"
#include "DDCheckerboard.h"
#include "RHSSrcPts.h"
#include "OutputVTK.h"
#include "GlobalIterativeSolverGMRES.h"
#include "FunctionCSV.h"
#include "FunctionConstant.h"
#include "Assembler.h"
#include "AssemblerHelmholtz2DFD.h"
#include "LocalLinearSolver.h"
#include "LocalLinearSolverEigenSparseLU.h"
#include "LocalLinearSolverPardiso.h"
#include "PML.h"
#include "PMLCubic2D.h"


using namespace std;
using namespace libconfig;

// The standard constructor
Input::Input(): mFileOpened(false)
{

}

// The copy constructor
Input::Input( const Input &other ): mFileOpened( false )
{
  cerr << "WARNING: Input file not opened while assigning!" << endl;
}

// The destructor
Input::~Input()
{

}

// The assignment operator
Input &Input::operator=( const Input &other )
{
  // Only assign if other is dfferent from this object
  if( this != &other )
  {
    cerr << "WARNING: Input file not opened while assigning!" << endl;
    mFileOpened = false;
  }

  return *this;
}

// The funcitont hat opens the file
void Input::openFile( string name )
{
  // Read the file. If there is an error, report it and exit.
  try
  {
    mConfig.readFile(name.c_str());
  }
  catch(const FileIOException &fioex)
  {
    std::cerr << "I/O error while reading file." << std::endl;
    exit(1);
  }
  catch(const ParseException &pex)
  {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
      << " - " << pex.getError() << std::endl;
    exit(1);
  }

  mFileOpened = true;
}

// The function returning the number of problems
UInt Input::getNumProblems() const
{
  // Make sure that there is a file opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // Get the root handle and make sure that problems are defined in the file
  const Setting &root = mConfig.getRoot();
  if( !root.exists("Problems") )
  {
    cerr << "ERROR: No Problems defined in input file" << endl;
    exit(1);
  }

  // return the number of problems
  return root["Problems"].getLength();
}

// The funciton reading in a problem
void Input::readProblemDD( DD *&dd, UInt probID ) const
{
  // Make sure there is a file opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // Get the root handle and make sure there are problems defined
  const Setting &root = mConfig.getRoot();
  if( !root.exists("Problems") )
  {
    cerr << "ERROR: No Problems defined in input file" << endl;
    exit(1);
  }

  // make sure that the problem ID exists
  if( probID >= root["Problems"].getLength() )
  {
    cerr << "ERROR: Invalid problem ID" << endl;
    exit(1);
  }

  // get the problem handle
  const Setting &problem = root["Problems"][probID];
  // get the DDType
  string DDType;
  if( !problem.lookupValue("DomainDecomposition", DDType ) )
  {
    cerr << "ERROR: No DomainDecomposition prescribed in problem!" << endl;
    exit(1);
  }
  // get the ID of the DD
  UInt DDID;
  if( !problem.lookupValue("DomainDecompositionID", DDID ) )
  {
    cerr << "ERROR: No DomainDecompositionID prescribed in problem!" << endl;
    exit(1);
  }

  // Make sure that the DDType is valid and read in the DD
  if( DDType=="CheckerboardDD" )
  {
    readCheckerboardDD( DDID, dd );
  }
  else
  {
    cerr << "ERROR: Invalid DomainDecomposition prescribed" << endl;
    exit(1);
  }
}

// The function that reads in the RHS
void Input::readProblemRHS( const DD *dd, RHS *&rhs, UInt probID ) const
{
  // Make sure that a file is opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // Get the root handle and make sure problems exist
  const Setting &root = mConfig.getRoot();
  if( !root.exists("Problems") )
  {
    cerr << "ERROR: No Problems defined in input file" << endl;
    exit(1);
  }

  // Make sure that the problem ID is valid
  if( probID >= root["Problems"].getLength() )
  {
    cerr << "ERROR: Invalid problem ID" << endl;
    exit(1);
  }

  // get the problem handle
  const Setting &problem = root["Problems"][probID];
  // get the RHSType
  string RHSType;
  if( !problem.lookupValue("RHS", RHSType ) )
  {
    cerr << "ERROR: No RHS prescribed in problem!" << endl;
    exit(1);
  }

  // make sure that there are RHSs defined
  if( !root.exists("RHSs") )
  {
    cerr << "ERROR: No RHSs prescribe in input file!" << endl;
    exit(1);
  }
  if( !root["RHSs"].exists( RHSType.c_str() ) )
  {
    cerr << "ERROR: Invalid RHS type prescribed" << endl;
    exit(1);
  }

  // get the handle on the sources
  const Setting &srcs = root["RHSs"][RHSType.c_str()];
  // allocate memory for the point sources
  vector<DoubleVector> srcPts(srcs.getLength());
  // loop over all sources and read them in
  for( UInt s=0; s<srcs.getLength(); s++ )
  {
    // get the dimension
    UInt dim;
    if( !srcs[s].lookupValue("Dim",dim) )
    {
      cerr << "ERROR: No Dim prescribed for source point" << endl;
      exit(1);
    }

    // get the x-value
    Double x;
    if( !srcs[s].lookupValue("x",x) )
    {
      cerr << "ERROR: No x prescribed in source point" << endl;
      exit(1);
    }

    // get the y-value
    Double y;
    if( !srcs[s].lookupValue("y",y) )
    {
      cerr << "ERROR: No y prescribed in source point" << endl;
      exit(1);
    }

    // set up the source point
    srcPts[s] = DoubleVector(dim);
    srcPts[s][0] = x;
    srcPts[s][1] = y;
  }

  // set up the RHS based on source points
  RHSSrcPts rhsSrcPts;
  // set source points
  rhsSrcPts.setSrcPts( srcPts );
  // set the epsilon used to approximate the point sources by
  // exponential functions to the mesh size
  rhsSrcPts.setEpsilon( dd->getMeshSize() ); 

  // copy the set up RHS into the result
  rhs = rhsSrcPts.copy();
}

// The function that reads in a GlobalIterativeSolver
void Input::readProblemSolver( GlobalIterativeSolver *&solver, UInt probID ) const
{
  // Make sure a file is opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // Get the root handle and make sure problems are defined
  const Setting &root = mConfig.getRoot();
  if( !root.exists("Problems") )
  {
    cerr << "ERROR: No Problems defined in input file" << endl;
    exit(1);
  }

  // Make sure that the problem ID is valid
  if( probID >= root["Problems"].getLength() )
  {
    cerr << "ERROR: Invalid problem ID" << endl;
    exit(1);
  }

  // get the handle to the problem
  const Setting &problem = root["Problems"][probID];
  // get the solver type
  string solverType;
  if( !problem.lookupValue("Solver", solverType ) )
  {
    cerr << "ERROR: No Solver prescribed in problem!" << endl;
    exit(1);
  }

  // read in the solver corresponding to the solverType
  readSolver( solverType, solver );
}

// The function reading in the Output
void Input::readProblemOutput( Output *&output, UInt probID ) const
{
  // Make sure a file is opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // Get the root handle and make sure problems are defined
  const Setting &root = mConfig.getRoot();
  if( !root.exists("Problems") )
  {
    cerr << "ERROR: No Problems defined in input file" << endl;
    exit(1);
  }

  // Make sure that the problem ID is valid
  if( probID >= root["Problems"].getLength() )
  {
    cerr << "ERROR: Invalid problem ID" << endl;
    exit(1);
  }

  // get the handle to the problem
  const Setting &problem = root["Problems"][probID];
  // get the output type
  string outputType;
  if( !problem.lookupValue("Output", outputType ) )
  {
    cerr << "ERROR: No Output prescribed in problem!" << endl;
    exit(1);
  }

  // use the outputType to read in the Output
  readOutput( outputType, output );
}

// The function reading in a particular GlobalIterativeSolver 
void Input::readSolver( string solverType, GlobalIterativeSolver *&solver ) const
{
  // Make sure a file is opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // Get the root handle and make sure GlobalSolvers are defined
  const Setting &root = mConfig.getRoot();
  if( !root.exists("GlobalSolvers") )
  {
    cerr << "ERROR: No GlobalSolvers defined in input file" << endl;
    exit(1);
  }

  // Make sure that the given solverType exists
  if( !root["GlobalSolvers"].exists(solverType.c_str()) )
  {
    cerr << "ERROR: Invalid Solver type prescribed!" << endl;
    exit(1);
  }

  // get the handle to the solver
  const Setting &sol = root["GlobalSolvers"][solverType.c_str()];
  // get the type of the solver
  string type;
  if( !sol.lookupValue("Type",type) )
  {
    cerr << "ERROR: No Type prescribed in Solver" << endl;
    exit(1);
  }
  // make sure that the type is valid
  // if so, read in the solver
  if( type=="GMRES" )
  {
    // read in the GMRES tolerance
    // if none is prescribed, use the standard value of 1e-6
    Double tol=1.0e-6;
    sol.lookupValue("Tol",tol); 

    // read in the maximum number of iterations
    // if maxIter=0 the GMRES method converges to the prescribed tolerance
    UInt maxIter=0;
    sol.lookupValue("MaxIter",maxIter);

    // read in if the GMRES method should be verbose
    // if verbose is true, then information is printed out as the system is solved
    bool verbose=true;
    sol.lookupValue("Verbose",verbose);

    // set up the GlobalIterativeSolver fo the GMRES method
    GlobalIterativeSolverGMRES temp;

    // set the read in values
    if( tol!=0.0 )
      temp.setTol( tol );
    if( maxIter!=0 )
      temp.setMaxIter( maxIter );
    if( verbose!=true )
      temp.setVerbose( verbose );

    // copy the object into the result
    solver = temp.copy();
  }
  else
  {
    cerr << "ERROR: Unknown Solver type!" << endl;
    exit(1);
  }
}

// The function reading in a particular output
void Input::readOutput( string outputType, Output *&output ) const
{
  // if the output type is empty, nothing should be done and we set the Output to NULL
  if( outputType=="" )
  {
    output = NULL;
    return;
  }

  // Make sure a file is opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // get the root handle and make sure Outputs are defined
  const Setting &root = mConfig.getRoot();
  if( !root.exists("Outputs") )
  {
    cerr << "ERROR: No Outputs defined in input file" << endl;
    exit(1);
  }
  // Make sure the output type exists
  if( !root["Outputs"].exists(outputType.c_str()) )
  {
    cerr << "ERROR: Invalid output type prescribed!" << endl;
    exit(1);
  }

  // get the handle to the output
  const Setting &out = root["Outputs"][outputType.c_str()];
  // read the output type
  string type;
  if( !out.lookupValue("Type",type ) )
  {
    cerr << "ERROR: No Type prescribed for Output" << endl;
    exit(1);
  }

  // make sure that the output type exists
  // if so, read in the output
  if( type=="VTK" )
  {
    // get the filename that should be used, the filename does not include the .vtk ending
    string file;
    if( !out.lookupValue("Filename",file) )
    {
      cerr << "ERROR: No Filename prescribed!" << endl;
      exit(1);
    }

    // set up the Output
    OutputVTK temp;
    // set the filename
    temp.setFilename(file);

    // copy the object into the result
    output = temp.copy();
  }
  else
  {
    cerr << "ERROR: Unknown output type!" << endl;
    exit(1);
  }
}

// The function read in a CheckerboardDD
void Input::readCheckerboardDD( UInt DDID, DD *&dd ) const
{
  // Make sure a file is opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // get the root handle
  const Setting &root = mConfig.getRoot();
  // Make sure DomainDecompositions are defined
  if( !root.exists("DomainDecompositions") )
  {
    cerr << "ERROR: No DomainDecompositions defined in input file" << endl;
    exit(1);
  }

  // get the handle to the DomainDecompositions
  const Setting &DomainDecompositions = root["DomainDecompositions"];
  if(!DomainDecompositions.exists("CheckerboardDD"))
  {
    cerr << "ERROR: No CheckerboardDDs defined in input file" << endl;
    exit(1);
  }

  // make sure the ID of the CheckerboardDD is valid
  if(DDID>= DomainDecompositions["CheckerboardDD"].getLength())
  {
    cerr << "ERROR: Invalid DomainDecompositionID prescribed!" << endl;
    exit(1);
  }

  // get the handle to the CheckerboardDD
  const Setting &CDD = DomainDecompositions["CheckerboardDD"][DDID];

  // get the dimension
  UInt dim;
  if(!CDD.lookupValue("Dim", dim))
  {
    cerr << "ERROR: No Dim prescribed in CheckerboardDD" << endl;
    exit(1);
  }
  if( dim != 2 )
  {
    cerr << "ERROR: Dim can only be 2!" << endl;
    exit(1);
  }
  
  // read in the domain
  string domain;
  if(!CDD.lookupValue("Domain",domain))
  {
    cerr << "ERROR: No Domain prescribed in CheckerboardDD" << endl;
    exit(1);
  }
  vector<DoubleVector> dom;
  readDomain(domain,dim,dom);

  // read in the assembler
  Assembler* assembler = NULL;
  // get the pde type
  string pde;
  if( !CDD.lookupValue("PDE", pde ) )
  {
    cerr << "ERROR: No PDE prescribed in CheckerboardDD" << endl;
    exit(1);
  }

  // if the pde type is valid, read in the corresponding assembler
  if( pde=="Helmholtz" )
  {
    if( dim==2 )
      assembler = readAssemblerHelmholtz2DFD(CDD,dom);
    else
    {
      cerr << "ERROR: Dim can only be 2!" << endl;
      exit(1);
    }
  }
  else
  {
    cerr << "ERROR: Unknown PDE prescribed!" << endl;
    exit(1);
  }

  // read in the number of elements used in each cell 
  UInt nX, nY;
  if( CDD.exists("NX") && CDD.exists("NY") )
  {
    CDD.lookupValue("NX",nX); 
    CDD.lookupValue("NY",nY);
  }
  else
  {
    cerr << "ERROR: Invalid way to prescribe N's" << endl;
    exit(1);
  }

  // read in the number of cells in the x direction
  UInt nCellsX;
  if( !CDD.lookupValue("NCellsX",nCellsX) )
  {
    cerr << "ERROR: No NCellsX prescribed in CheckerboardDD" << endl;
    exit(1);
  }
  // read in the number of cells in the y direction
  UInt nCellsY;
  if( !CDD.lookupValue("NCellsY",nCellsY) )
  {
    cerr << "ERROR: No NCellsY prescribed in CheckerboardDD" << endl;
    exit(1);
  }

  // define the local solver
  LocalLinearSolver *solver;
  // get the solver type
  string solverType;
  if( !CDD.lookupValue("Solver",solverType) )
  {
    cerr << "ERROR: No solver prescribed!" << endl;
    exit(1);
  }
  // set up the solvers if the solver type is known
  if( solverType=="Eigen" )
    solver = new LocalLinearSolverEigenSparseLU();
  else if( solverType=="Pardiso" )
    solver = new LocalLinearSolverPardiso();
  else
  {
    cerr << "ERROR: Unknown solver prescribed!" << endl;
    exit(1);
  }

  // set up the Checkerboard DD
  UIntVector nCells(dim);
  nCells[0] = nCellsX;
  nCells[1] = nCellsY;

  UIntVector nElsInCell(dim);
  nElsInCell[0] = nX;
  nElsInCell[1] = nY;

  dd = new DDCheckerboard();
  if( dim==2 )
    ((DDCheckerboard*)dd)->setUp( assembler, solver, dom, nCells, nElsInCell, MPI_COMM_WORLD );
  else
  {
    cerr << "ERROR: Only dim=2 is allowed!" << endl;
    exit(1);
  }

  // delete the temporary objects
  delete assembler;
  delete solver;
}

// The function reading in a domain
void Input::readDomain( string domainType, UInt dim, vector<DoubleVector> &dom ) const
{
  // Make sure a file is opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // get the root handle and make sure that domains are defined in the file
  const Setting &root = mConfig.getRoot();
  if( !root.exists("Domains") )
  {
    cerr << "ERROR: No Domains defined in input file" << endl;
    exit(1);
  }
  // get the handle on the domain and make sure that the domainType is found
  const Setting &domains = root["Domains"];
  if( !domains.exists( domainType.c_str() ) )
  {
    cerr << "ERROR: Unknown Domain prescribed" << endl;
    exit(1);
  }
  const Setting &domain = domains[domainType.c_str()];

  // get the dimension
  UInt domainDim;
  if( !domain.lookupValue( "Dim", domainDim ) )
  {
    cerr << "ERROR: No Dim prescribed in Domain" << endl;
    exit(1);
  }

  // make sure that the dimension is consistent
  if( dim!=domainDim )
  {
    cerr << "ERROR: Inconsistent dimensions in Domain" << endl;
    exit(1);
  }

  // get the domain interval in the x-direction
  // left boundary
  Double ax;
  if( !domain.lookupValue( "ax", ax ) )
  {
    cerr << "ERROR: No ax prescribed for " << domainType << endl;
    exit(1);
  }
  // right boundary
  Double bx;
  if( !domain.lookupValue( "bx", bx ) )
  {
    cerr << "ERROR: No bx prescribed for " << domainType << endl;
    exit(1);
  }

  // get the domain interval in the y-direction
  // bottom boundary
  Double ay;
  if( !domain.lookupValue( "ay", ay ) )
  {
    cerr << "ERROR: No ay prescribed for " << domainType << endl;
    exit(1);
  }
  // top boundary
  Double by;
  if( !domain.lookupValue( "by", by ) )
  {
    cerr << "ERROR: No by prescribed for " << domainType << endl;
    exit(1);
  }

  // set up the domain from the read values
  dom = vector<DoubleVector>(dim);
  dom[0] = DoubleVector(2);
  dom[0][0] = ax;
  dom[0][1] = bx;
  dom[1] = DoubleVector(2);
  dom[1][0] = ay;
  dom[1][1] = by;
}

// The function that reads in a squared slowness
Function *Input::readSquaredSlowness( string ssName ) const
{
  // Make sure a file is opened
  if( !mFileOpened )
  {
    cerr << "ERROR: No file opened for input" << endl;
    exit(1);
  }

  // get the root handle
  const Setting &root = mConfig.getRoot();
  // make sure that SquaredSlownesses are defined
  if( !root.exists("SquaredSlownesses") )
  {
    cerr << "ERROR: No SquaredSlownesses defined in input file" << endl;
    exit(1);
  }
  // make sure that the squared slowness exists
  if( !root["SquaredSlownesses"].exists(ssName.c_str()) )
  {
    cerr << "ERROR: Invalid SquaredSlowness prescribed" << endl;
    exit(1);
  }

  // get the handle on the squared slowness
  const Setting &ss = root["SquaredSlownesses"][ssName.c_str()];

  // get the type
  string type;
  if( !ss.lookupValue("Type",type) )
  {
    cerr << "ERROR: No Type prescribed for the squared slowness" << endl;
    exit(1);
  }
  // if the type is valid, read in the appropriate squared slowness
  if( type == "Constant" )
  {
    // constant squared slowness
    FunctionConstant fun;
    // we define it to be 1 so that the effective frequency is equal to omega
    fun.setVal(1.0);
    return fun.copy();  
  }
  else if( type == "CSV" )
  {
    // get the filename
    string file;
    if( !ss.lookupValue("Filename",file) )
    {
      cerr << "ERROR: No filename prescribed in CSV squared slowness" << endl;
      exit(1);
    }
    // get the dimension
    UInt dim;
    if( !ss.lookupValue("Dim",dim) )
    {
      cerr << "ERROR: No Dim prescribed in CSV squared solwness" << endl;
      exit(1);
    }
    // get the domain of the matrix defined in the file
    string domain;
    if( !ss.lookupValue("Domain",domain) )
    {
      cerr << "ERROR: No Domain prescribed in CSV squared slowness" << endl;
      exit(1);
    }
    // read the prescribed domain
    vector<DoubleVector> dom;
    readDomain(domain,dim,dom);

    // set up the CSV function
    FunctionCSV fun;
    fun.setFile( file, dom );
    return fun.copy();
  }
  else
  {
    cerr << "ERROR: Invalid SquaredSlowness type prescribed" << endl;
    exit(1);
  }

  // if there was a mistake, return NULL
  return NULL;
}

// The function that reads in a 2D finite difference Helmholtz assembler
Assembler *Input::readAssemblerHelmholtz2DFD(const Setting &CDD, const vector<DoubleVector> &domain ) const
{
  // define the assembler
  AssemblerHelmholtz2DFD assembler;

  // read the frequency omega
  Double omega;
  if(!CDD.lookupValue("Omega",omega))
  {
    cerr << "ERROR: No omega prescribed!" << endl;
    exit(1);
  }
  cout << "omega=" << omega << endl;
  // read the PML width given by the minimum number of wavelengths
  Double minWavelengthsInPML;
  if(!CDD.lookupValue("MinWavelengthsInPML",minWavelengthsInPML))
  {
    cerr << "ERROR: No MinWavelengthsInPML prescribed!" << endl;
    exit(1);
  }

  // get the squared slowness type that should be used in the assembler
  string ssType;
  if( !CDD.lookupValue("SquaredSlowness",ssType) )
  {
    cerr << "ERROR: No SquaredSlowness prescribed!" << endl;
    exit(1);
  }
  // set up the squred slowness
  Function *m = readSquaredSlowness( ssType );

  // get the discretization order
  UInt order;
  if( !CDD.lookupValue("Order",order) )
  {
    cerr << "ERROR: No Order prescribed in CheckerboardDD" << endl;
    exit(1);
  }

  // define the PML
  PMLCubic2D pml;
  // set the absorption constant
  Double C = log(omega);
  pml.setAbsorptionCoeff(C);
  cout << C << endl;
  // determine the wavelength and compute the maximum wavelength for the given squared slowness and frequency
  Double wavelength1 = (domain[0][1]-domain[0][0])*2*pi/omega;
  Double wavelength2 = (domain[1][1]-domain[1][0])*2*pi/omega;
  Double wavelength = min(wavelength1,wavelength2);
  Double minimum = m->getMinReal();
  wavelength /= sqrt(minimum);
  cout << wavelength << endl;

  // update the PML width
  Double pmlWidth = minWavelengthsInPML*wavelength;
  pml.setWidth(pmlWidth);

  // use the read and adjusted data to set up the assembler
  assembler.setPML(&pml);
  assembler.setOmega(omega);
  m->setDomain(domain);
  assembler.setM(m);
  assembler.setOrder(order);
  // return a pointer to a copy of the set up assembler
  return assembler.copy();
}
