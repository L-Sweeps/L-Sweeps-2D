#ifndef PROBLEM_H
#define PROBLEM_H

#include "typedef.h"


class DD;
class RHS;
class GlobalIterativeSolver;
class Output;
class Input;

//! The class that collects all required objects for a problem and uses them to solve the problem
class Problem
{
  protected:
    //! A pointer to the DD class
    DD *mDD;
    //! A pointer to the RHS class
    RHS *mRHS;
    //! A pointer to the GlobalSolver class
    GlobalIterativeSolver *mSolver;
    //! A pointer to the Output class
    Output *mOutput;

  public:
    //! The standard constructor, constructs an empty problem class
    Problem();
    
    //! The copy constructor
    Problem( const Problem &other );
    
    //! The destructor, frees all memory
    ~Problem();

    //! The assignment operator
    Problem &operator=( const Problem &other );

    //! The function that reads a problem from an Input class
    /*!
     * \p input The Input class
     * \p probID The ID of the problem in the input file that should be read
     */
    void read( const Input &input, UInt probID );

    //! The function that solves the problem
    void solve() const;
};

#endif
