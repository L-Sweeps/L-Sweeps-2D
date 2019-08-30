#ifndef LOCALLINEARSOLVEREIGENSPARSELU_H
#define LOCALLINEARSOLVEREIGENSPARSELU_H

#include "typedef.h"
#include "LocalLinearSolver.h"
#include <Eigen/SparseLU>

using namespace::Eigen;

//! A class that realizes a local solver using a direct solver from Eigen
class LocalLinearSolverEigenSparseLU : public LocalLinearSolver
{
  protected:
    //! The system matrix of the system to be solved
    const ComplexSparseMatrix *mA;

    //! The solver object, this is just the eigen solver object
    SparseLU<ComplexSparseMatrix> mSolver;

    //! The flag if the system is factorized
    bool mIsFactorized;

  public:
    //! The standard constructor, generates an empty object
    LocalLinearSolverEigenSparseLU(): mA(NULL), mSolver(), mIsFactorized(false) {}
    
    //! The destructor
    virtual ~LocalLinearSolverEigenSparseLU() {}

    //! The copy function
    /*!
     * Generates a new copy of this object and returns a pointer to the copy.
     */
    virtual LocalLinearSolver *copy() const;

    //! The function that sets the matrix of the system to be solved
    /*!
     * \param A The matrix to be set, it has to be given in compressed form
     */
    virtual void setMatrix( ComplexSparseMatrix *A );
    
    //! The function to factorize the system
    virtual void factorize();
    
    //! The function to solve the system, solves the system using Eigen
    /*!
     * \param b The right hand side vector
     * \param res The solution vector
     */
    virtual void solve( ComplexVector &b, ComplexVector &res ) const;

  private:
    //! The (disabled) copy constructor, should not be called from the outside
    /*!
     * The fact that it shouldn't be called from the outside is just inherited
     * from the Eigen library which has the copy constructor of the solver object
     * disabled.
     */
    LocalLinearSolverEigenSparseLU(const LocalLinearSolverEigenSparseLU &other);

    //! The (disabled) assignment, should not be called from the outside
    /*!
     * The fact that it shouldn't be called from the outside is just inherited
     * from the Eigen library which has the assignment operator of the solver object
     * disabled.
     */
    LocalLinearSolverEigenSparseLU &operator=(const LocalLinearSolverEigenSparseLU &other);

};

#endif
