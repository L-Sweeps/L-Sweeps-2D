#ifndef GLOBALITERATIVESOLVER_H
#define GLOBALITERATIVESOLVER_H

#include "DD.h"
#include <vector>

using namespace std;

//! The abstract class for a global iterative solver 
/*!
 * This class is overloaded by specific solver routines. This is a solver
 * that uses the domain decomposition preconditioner and system matrix information
 * and uses them in an iterative solver in order to solve the global system defined
 * in the domain decomposition.
 */
class GlobalIterativeSolver
{
  protected:
    //! The parent domain decomposition
    const DD *mParent;

  public:
    //! The standard constructor
    GlobalIterativeSolver();

    //! The copy constructor
    GlobalIterativeSolver( const GlobalIterativeSolver &other );

    //! The destructor
    virtual ~GlobalIterativeSolver() {};

    //! The assignment operator
    GlobalIterativeSolver &operator=( const GlobalIterativeSolver &other );
    
    //! The function that sets the parent domain decomposition
    void setParent( const DD *parent );

    //! The copy function
    /*!
     * This function copies the (derived) object properly and returns a pointer to
     * the abstract class. This is important to be able to store different types
     * of assembly routines in other parts of the code while not changing the code
     * and consequently makes the code more modular.
     */
    virtual GlobalIterativeSolver* copy() const = 0;

    //! The function that computes a global solution
    virtual void solve(const vector<ComplexVector> &rhs, vector<ComplexVector> &res) const = 0;
};

#endif
