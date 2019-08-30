#ifndef DDCHECKERBOARD_H
#define DDCHECKERBOARD_H

#include "DD.h"
#include "typedef.h"
#include "Subdomain.h"
#include <mpi.h>
#include "Assembler.h"


class Function;

//! A class that implemets a checkerboard domain decompositon (DD)
class DDCheckerboard : public DD
{
  protected:
    //! The local cells in each process
    vector<Subdomain*> mCells;
    
    //! The matrx that saves which subdomain is processed by what process
    UIntMatrix mIndToProc;

    //! The indices determining the interior DOFs
    /*!
     * The indices are given as a vector of vectors. The outer vector has lenght equal
     * to the number of cells, the inner vector holds the indices corresponding to the
     * subdomain. The indices are based on all DOFs of the local problem and hold the indices
     * corresponding to the interior DOFs (\fs \boldsymbol{\Omega}_{ij}\fs in the paper)
     */
    vector<UIntVector> mValsToKeep;

    //! the indices determining the interior DOFs with shifted traces in both direction.
    /*!
     * \sa mValsToKeep
     */
    vector<UIntVector> mValsToKeepShifted;
    
    //! the indices determining the interior DOFs with shifted traces in the Left-Right direction.
    /*!
     * \sa mValsToKeep
     */
    vector<UIntVector> mValsToKeepShiftedLR;
    
    //! the indices determining the interior DOFs with shifted traces in the Up-Down direction.
    /*!
     * \sa mValsToKeep
     */
    vector<UIntVector> mValsToKeepShiftedUD;

    //! The number of cells in the row direction
    UInt mNCells1;
    
    //! The number of cells in the column direction
    UInt mNCells2;

  public:
    //! The standard constructor
    DDCheckerboard() {}
    
    //! The copy constructor
    DDCheckerboard( const DDCheckerboard &other );
    
    //! The destructor
    ~DDCheckerboard();

    //! The copy function
    virtual DD* copy() const;

    //! The assignment operator
    DDCheckerboard &operator=( const DDCheckerboard &other );

    //! The function that sets up a checkerboard domain decomposition for a given assmbler and solver 
    /*!
     * \param assembler The assembler used to discretize the local problems \sa Assembler.
     * \param solver The direct solver used to solve the local linear systems \sa LocalLinearSolver.
     * \param domain The domain of interest (without the PMLs), stored in each direction, e.g.
     *               the unit square is saved as [ [0,1], [0,1] ]
     * \param nCells The number of cells in each direction
     * \param nElsInCell The number of elements that should be used inside (between trace1 and traceN) in each direction
     * \param comm The MPI communicator
     */
    void setUp( const Assembler *assembler, const LocalLinearSolver* solver, vector<DoubleVector> domain, UIntVector nCells, UIntVector nElsInCell, MPI_Comm comm );

    //! The function that factorizes each subdomain
    void factorize();

    //! The function that assembles the rhs in each subdomain for a given RHS 
    /*!
     * \param rhs The RHS to be used
     * \param res The vector that stores the result, is overwritten
     */
    void assembleRHS( const RHS *rhs, vector<ComplexVector> &res ) const;
    
    //! The function that applies the preconditioner (L-Sweeps)
    /*!
     * \param x The vector the preconditioner should be applied to
     * \param res The vector that stores the result
     */
    void applyPrecond( const vector<ComplexVector> &x, vector<ComplexVector> &res ) const;

    //! The function that appies the local system matrix to a given vector
    /*!
     * \param x The vector the system matrix should be applied to.
     *          This is stored as a vector of vectors in each subdomain. Each subdomain 
     *          only stores its local indices (the ones that are given in mValsToKeep, 
     *          for example not trace0 or traceNP)
     * \param The resulting vector, stored in the same way as \p x
     */
    void applySystemMatrix( const vector<ComplexVector> &x, vector<ComplexVector> &res ) const;

    //! The function returning the dimension of the underlying problem
    UInt getDim() const;

    //! The function returning total number of subdomains (cells)
    UInt getNCells() const;

    //! The function that returns a pointer to a cell (subdomain)
    /*!
     * \param i The index of the subdomain that should be returned
     */
    Subdomain* getCell(UInt i) const;
  private:
    //! The function that applies all L-Sweeps in all directions
    /*!
     * \param x The vector the sweep should be applied to. 
     *          This is stored as a vector of vectors in each subdomain. Each subdomain 
     *          only stores its local indices (the ones that are given in mValsToKeep, 
     *          for example not trace0 or traceNP)
     * \param The resulting vector, stored in the same way as \p x
     * \param shifted The flag to determine if the shifted or non-shifted sweeps should be applied
     */
    void applySweep( const vector<ComplexVector> &x, vector<ComplexVector> &res, bool shiftedLR, bool shiftedUD ) const;

    //! The function that returns the global mesh size
    Double getMeshSize() const;

    //! The function that takes a vector of values defined over all DOFs in each subdomain and reduces it to the subdomains that should be used from each subdomain
    /*!
     * \param x The vectors to be reduced
     * \param res The resulting vectors
     */
    void reduce( const vector<ComplexVector> &x, vector<ComplexVector> &res, bool shiftedLR = false, bool shiftedUD=false ) const;

    //! The function that takes a vector of reduced values (the values each subdomain should store) and expands it to all DOFs of each subdomain by zero extension
    /*!
     * \param x The vectors to be expanded
     * \param res The resulting vectors
     */
    void expand( const vector<ComplexVector> &x, vector<ComplexVector> &res, bool shiftedLR = false, bool shiftedUD=false ) const;

    //! The function that adds the results of the local solutions in each subdomain
    /*!
     * \param x The input vectors of the local source vectors, given as a vector of values corresponding to the interior DOFs (\fs \boldsymbole{\Omega}_{ij}\fs in the paper)
     * \param res The input vectors to which the contribution of the local solutions are added, given as a vector of values corresponding to the interior DOFs (\fs \boldsymbole{\Omega}_{ij}\fs in the paper)
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     * */
    void addLocSol( const vector<ComplexVector> &x, vector<ComplexVector> &res, vector<vector<ComplexVector>> &tracesLR, vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const;

    //! The function that adds the results of a right sweep to given vectors defined over the DD
    /*!
     * \param res The vector to which the contribution of the sweep should be added 
     * \param tracesLR The Left-Right traces that are extended to the right 
     * \param tracesUD The vector where the UpDown trace information is stored
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     */
    void addRightSweep( vector<ComplexVector> &res, const vector<vector<ComplexVector>> &tracesLR, vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const;

    //! The function that adds the results of a left sweep to given vectors defined over the DD
    /*!
     * \param res The vector to which the contribution of the sweep should be added 
     * \param tracesLR The Left-Right traces that are extended to the left 
     * \param tracesUD The vector where the UpDown trace information is stored
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     */
    void addLeftSweep( vector<ComplexVector> &res, const vector<vector<ComplexVector>> &tracesLR, vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const;


    //! The function that adds the results of an up sweep to given vectors defined over the DD
    /*!
     * \param res The vector to which the contribution of the sweep should be added 
     * \param tracesLR The vector where the Left-Right trace information is stored
     * \param tracesUD The Up-Down traces that are extended to the top
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     */
    void addUpSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> &tracesLR, const vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const;

    //! The function that adds the results of a down sweep to given vectors defined over the DD
    /*!
     * \param res The vector to which the contribution of the sweep should be added 
     * \param tracesLR The vector where the Left-Right trace information is stored
     * \param tracesUD The Up-Down traces that are extended to the bottom
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     */
    void addDownSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> &tracesLR, const vector<vector<ComplexVector>> &tracesUD, bool shiftedLR, bool shiftedUD ) const;

    //! The function that adds the results of an L-sweep from bottom left to top right to given vectors defined over the DD
    /*!
     * \param res The vector to which the contribution of the sweep should be added 
     * \param tracesLR The input LeftRight traces for the sweep, is changed in the function
     * \param tracesUD The input UpDown traces for the sweep, is changed in the function
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     */
    void addBL2TRSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> tracesLR, vector<vector<ComplexVector>> tracesUD, bool shiftedLR, bool shiftedUD ) const;
    
    //! The function that adds the results of an L-sweep from top right to bottom left to given vectors defined over the DD
    /*!
     * \param res The vector to which the contribution of the sweep should be added 
     * \param tracesLR The input LeftRight traces for the sweep, is changed in the function
     * \param tracesUD The input UpDown traces for the sweep, is changed in the function
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     */
    void addTR2BLSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> tracesLR, vector<vector<ComplexVector>> tracesUD, bool shiftedLR, bool shiftedUD ) const;
    
    //! The function that adds the results of an L-sweep from bottom right to top left to given vectors defined over the DD
    /*!
     * \param res The vector to which the contribution of the sweep should be added 
     * \param tracesLR The input LeftRight traces for the sweep, is changed in the function
     * \param tracesUD The input UpDown traces for the sweep, is changed in the function
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     */
    void addTL2BRSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> tracesLR, vector<vector<ComplexVector>> tracesUD, bool shiftedLR, bool shiftedUD ) const;
    
    //! The function that adds the results of an L-sweep from top left to bottom right to given vectors defined over the DD
    /*!
     * \param x The input vectors to which the L-sweep should be added
     * \param res The resulting vectors
     * \param tracesLR The input LeftRight traces for the sweep, is changed in the function
     * \param tracesUD The input UpDown traces for the sweep, is changed in the function
     * \param shifted The flag if the shifted or non-shifted traces should be used
     */
    void addBR2TLSweep( vector<ComplexVector> &res, vector<vector<ComplexVector>> tracesLR, vector<vector<ComplexVector>> tracesUD, bool shiftedLR, bool shiftedUD ) const;

    //! The function that communicates trace information from bottom left to top right
    /*!
     * \param xEx The input and output vector of values in each subdomain, the values are assumed to be given in expanded form \sa expand.
     * \param shifted The flag if the input/output vector holds the values in shifted form (or not)
     * This function takes the top/right trace (traceN) and communicates them to the top/right element.
     */
    void communicateBL2TRTraces( vector<ComplexVector> &xEx, bool shifted ) const;
    
    //! The function that communicates trace information from top right to bottom left
    /*!
     * \param xEx The input and output vector of values in each subdomain, the values are assumed to be given in expanded form \sa expand.
     * \param shifted The flag if the input/output vector holds the values in shifted form (or not)
     *
     * This function takes the bottom/left trace (traceN) and communicates them to the bottom/left element.
     */
    void communicateTR2BLTraces( vector<ComplexVector> &xEx, bool shifted ) const;

    //! The function that converts DOFs from shifted to non-shfited form in the Left-Right direction
    /*!
     * \param x The input/output vector, given in reduced form \sa reduce
     * \param shiftedUD A flag if the Up-Down traces are shifted
     */
    void convertShiftedToNonShiftedLR( vector<ComplexVector> &x, bool shiftedUD ) const;
    
    //! The function that converts DOFs from shifted to non-shfited form in the Up-Down direction
    /*!
     * \param x The input/output vector, given in reduced form \sa reduce
     * \param shiftedLR A flag if the Left-Right traces are shifted
     */
    void convertShiftedToNonShiftedUD( vector<ComplexVector> &x, bool shiftedLR ) const;
    
    //! The function that converts DOFs from non-shifted to shfited form in the Left-Right direction
    /*!
     * \param x The input/output vector, given in reduced form \sa reduce
     * \param shiftedUD A flag if the Up-Down traces are shifted
     */
    void convertNonShiftedToShiftedLR( vector<ComplexVector> &x, bool shiftedUD ) const;
    
    //! The function that converts DOFs from non-shifted to shfited form in the Up-Down direction
    /*!
     * \param x The input/output vector, given in reduced form \sa reduce
     * \param shiftedLR A flag if the Left-Right traces are shifted
     */
    void convertNonShiftedToShiftedUD( vector<ComplexVector> &x, bool shiftedLR ) const;

    //! The function that separates the RHS into one of the four RHS windowed for each (shifted) DD
    /*!
     * \param rhs The vector to be separated
     * \param res The results vector
     * \param shiftedLR The flag determining if the traces should be shifted in the Left-Right direction
     * \param shiftedUD The flag determining if the traces should be shifted in the Up-Down direction
     */
    void getSeparatedRHS(const vector<ComplexVector> &rhs, vector<ComplexVector> &res, bool shiftedLR, bool shiftedUD ) const;
};

#endif
