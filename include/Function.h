#ifndef FUNCTION_H 
#define FUNCTION_H

#include "typedef.h"

//! The abstract class for a function
/*!
 * The class is overloaded by specific implementaitons of functions.
 * Examples where these classes are used are for the squared slowness or the sources.
 */
class Function
{
  public:
    //! The standard constructor
    Function();

    //! The copy constructor
    Function( const Function &other );

    //! The destructor
    virtual ~Function();

    //! The assignment operator
    Function &operator=( const Function &other );

    //! The copy function
    /*!
     * This function copies the (derived) object properly and returns a pointer to
     * the abstract class. This is important to be able to store different types
     * of assembly routines in other parts of the code while not changing the code
     * and consequently makes the code more modular.
     */
    virtual Function *copy() const = 0;

    //! The function that clears all the cache of the function
    virtual void clearCache() const = 0;

    //! The function that returns the minimum of the real part of the function
    virtual Double getMinReal() const = 0;
    
    //! The function that returns the maximum of the real part of the function
    virtual Double getMaxReal() const = 0;
    
    //! The function that returns the minimum of the imaginary part of the function
    virtual Double getMinImag() const = 0;
    
    //! The function that returns the maximum of the imaginary part of the function
    virtual Double getMaxImag() const = 0;

    //! The function that sets the domain of definition of the functoin
    /*!
     * The domain of definition is assumed to be the square.
     *
     * \p dom The domain
     */
    virtual void setDomain( const vector<DoubleVector> &dom ) const = 0;

    //! The function that returns the domain of definition
    /*!
     * \p res The result
     */
    virtual void getDomain( vector<DoubleVector> &res ) const = 0;

    //! The function that evaluates the Function
    /*!
     * \p x The point at which the function is evaluated
     */
    virtual Complex eval( DoubleVector x ) const = 0;
};

#endif
