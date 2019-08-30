#ifndef RHS_H
#define RHS_H

#include "typedef.h"

//! The abstract class for a right hand side
/*!
 * This class has to be overwritten with specific implementations.
 */
class RHS
{
  public:
    //! The standard constructor
    RHS();

    //! The copy constructor
    RHS( const RHS &other );

    //! The destructor
    virtual ~RHS();

    //! The assignment operator
    RHS &operator=( const RHS &other );

    //! The copy funciton to properly copy an object
    virtual RHS *copy() const = 0;

    //! The function that evaluates a right hand side
    virtual Complex eval( DoubleVector x ) const = 0; 
};

#endif
