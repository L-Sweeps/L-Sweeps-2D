#ifndef FUNCTIONCONSTANT_H
#define FUNCTIONCONSTANT_H

#include "Function.h"

//! A class that realizes a constant (complex) function
class FunctionConstant : public Function
{
  protected:
    //! The value of the constant function
    Complex mVal;

  private:
    //! The domai of definition
    mutable vector<DoubleVector> mDom;

  public:
    //! The standard constructor
    /*!
     * Constructs an empty element. The set functions
     * should be used to set the information.
     */
    FunctionConstant();

    //! The copy constructor
    /*!
     * Properly copies all member variables.
     */
    FunctionConstant( const FunctionConstant &other );
    
    //! The destructor
    virtual ~FunctionConstant();

    //! The assignment operator
    FunctionConstant &operator=( const FunctionConstant &other );

    //! The function that clears the cache, this function does nothing here
    virtual void clearCache() const;

    //! The function that returns the minimum value of the real part
    virtual Double getMinReal() const;

    //! The function that returns the maximum value of the real part
    virtual Double getMaxReal() const;

    //! The function that returns the minimum value of the imaginary part
    virtual Double getMinImag() const;

    //! The function that returns the maximum value of the imaginary part
    virtual Double getMaxImag() const;

    //! The function that sets the domain of definition
    /*!
     * \p dom The domain of definition
     */
    virtual void setDomain( const vector<DoubleVector> &dom ) const;
    
    //! the function that returns the domain of definition
    /*!
     * \p res The result
     */
    virtual void getDomain( vector<DoubleVector> &res ) const;

    //! The function that sets the value of the constant function
    /*!
     * \p val The value to be set
     */
    void setVal( Complex val );

    //! The copy function
    virtual Function *copy() const;

    //! The function that evaluates the function
    /*!
     * \p x The point at which the function is evaluated
     */
    virtual Complex eval( DoubleVector x ) const;
};

#endif
