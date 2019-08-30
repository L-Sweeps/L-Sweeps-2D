#ifndef INPUT_H
#define INPUT_H

#include "DD.h"
#include "Assembler.h"
#include "RHS.h"
#include "GlobalIterativeSolver.h"
#include "Output.h"
#include "Problem.h"
#include <string>
#include "libconfig.h++"
#include "Function.h"

using namespace std;
using namespace libconfig;

//! The class that handles the input file
class Input
{
  protected:
    //! The config object of libconfig
    Config mConfig;

    //! The flag determining if the file is opened
    bool mFileOpened;

  public:
    //! The standard constructor
    /*!
     * Sets up an empty object.
     */
    Input();

    //! The copy constructor
    Input( const Input &other );

    //! The destructor
    ~Input();

    //! The assignment operator
    Input &operator=( const Input &other );

    //! The function to open a file
    void openFile( string name );

    //! The function to read in a domain decomposition
    /*!
     * \p dd The pointer to the constructed DD, \sa DD
     * \p probID The ID of the problem from which the DD is read
     */
    void readProblemDD( DD *&dd, UInt probID ) const;
    
    //! The funciton that reads the RHS info
    /*!
     * \p dd The DD serving as the parent of the rhs, is input, \sa DD
     * \p rhs The RHS read in, \sa RHS
     * \p probID The ID of the problem from which the DD is read
     */
    void readProblemRHS( const DD* dd, RHS *&rhs, UInt probID ) const;
    
    //! The function that reads in solver information of the Global solver
    /*!
     * \p solver The pointer to the solver read in, \sa GlobalIterativeSolver
     * \p probID The ID of the problem from which the DD is read
     */
    void readProblemSolver( GlobalIterativeSolver *&solver, UInt probID ) const;
    
    //! The function that reads the output info
    /*!
     * \p output The pointer to the output read in, \sa Output
     * \p probID The ID of the problem from which the DD is read
     */
    void readProblemOutput( Output *&output, UInt probID ) const;

    //! The function that returns the number of problems defined in the file
    UInt getNumProblems() const;

  private:
    //! The function reading in a CheckerboardDD
    /*!
     * \p DDID The ID of the DD defined in the file that should be read in
     * \p dd The CheckerboardDD read in, \sa CheckerboardDD
     */
    void readCheckerboardDD( UInt DDID, DD *&dd ) const;

    //! The function reading in the AssemblerHelmholtz2DFD
    /*!
     * \p CDD The Setting handle used to read the assembler
     * \p domain The domain on which the assembler is defined
     * \return The pointer to the read in AssemblerHelmholtz2DFD class
     */
    Assembler *readAssemblerHelmholtz2DFD( const Setting &CDD, const vector<DoubleVector> &domain ) const;
    
    //! The function reading in a domain
    /*!
     * \p domainType The type of the domain
     * \p dim The dimension of the domain
     * \p dom The domain read in
     */
    void readDomain( string domainType, UInt dim, vector<DoubleVector> &dom ) const;
    
    //! The function reading in the squared slowness as a Function, \sa Function
    /*!
     * \p ssType The type of the squared slowness
     * \return The pointer to the read in Function
     */
    Function *readSquaredSlowness( string ssType ) const;
    
    //! The function reading in the GlobalIterativeSolver, \sa GlobalIterativeSolver
    /*!
     * \p solvertype The type of the solver
     * \p  The pointer to the read in solver
     */
    void readSolver( string solvertype, GlobalIterativeSolver *&solver ) const ;
    
    //! The function reading in the Output, \sa Output
    /*!
     * \p outputtype The type of the Output
     * \p The pointer to the read in output
     */
    void readOutput( string outputtype, Output *&output ) const ;
};

#endif
