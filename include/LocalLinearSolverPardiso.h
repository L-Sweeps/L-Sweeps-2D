#ifndef LOCALLINEARSOLVERPARDISO_H
#define LOCALLINEARSOLVERPARDISO_H

#include "typedef.h"
#include "LocalLinearSolver.h"

using namespace::Eigen;

//! A class that realizes a local solver using a direct solver from Eigen
class LocalLinearSolverPardiso : public LocalLinearSolver
{
  protected:
    //! The system matrix of the system to be solved
    ComplexSparseMatrix *mA;

    mutable void *mPt[64];
    mutable int mIParm[64];
    mutable Double mDParm[64];

    int *mIA;
    int *mJA;

    //! The flag if the system is factorized
    bool mIsFactorized;
    bool mHasMatrix;

  public:
    //! The standard constructor, generates an empty object
    LocalLinearSolverPardiso();
    
    //! The destructor
    virtual ~LocalLinearSolverPardiso();

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
    LocalLinearSolverPardiso(const LocalLinearSolverPardiso &other);

    //! The (disabled) assignment, should not be called from the outside
    /*!
     * The fact that it shouldn't be called from the outside is just inherited
     * from the Eigen library which has the assignment operator of the solver object
     * disabled.
     */
    LocalLinearSolverPardiso &operator=(const LocalLinearSolverPardiso &other);

    void clear();

};

#endif
