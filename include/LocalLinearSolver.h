#ifndef LOCALLINEARSOLVER_H
#define LOCALLINEARSOLVER_H

//! The abstract class for a local solver
/*!
 * This class is overloaded by specific routines for local solvers. The abstract class is used in other parts of the code (e.g. a domain decomposition).
 */
class LocalLinearSolver
{
  public:
    //! The standard constructor
    LocalLinearSolver() {}

    //! The copy constructor
    LocalLinearSolver( const LocalLinearSolver &other ) {}
    
    //! The destructor
    virtual ~LocalLinearSolver() {}

    //! The copy function
    /*!
     * This function copiesthe (derived) object properly and returns a pointer to
     * the abstract class. This is important to be able to store different types
     * of local linear solver routines in other parts of the code while not changing
     * the code and consequently makes the code more modular.
     */
    virtual LocalLinearSolver *copy() const = 0;

    //! A function setting the matrix of the local linear solver
    virtual void setMatrix( ComplexSparseMatrix *A ) = 0;

    //! A function that factorizes the matrix assigned to the solver
    virtual void factorize() = 0;

    //! A function that solves the stored system for a rhs
    /*!
     * \param b The right hand side vector
     * \param res The solution vector
     */
    virtual void solve( ComplexVector &b, ComplexVector &res ) const = 0;
};

#endif
