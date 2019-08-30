#ifndef RHSSRCPTS_H
#define RHSSRCPTS_H

#include "RHS.h"

//! The clas that evaluates a right hand side consisting of one or more point sources
/*!
 * This class approximates the point sources with exponential functions.
 */
class RHSSrcPts : public RHS
{
  protected:
    //! The locations of the point sources
    vector<DoubleVector> mSourcePts;

    //! The accuracy to which the point sources should be approximated by exponential functions
    Double mEps;

  public:
    //! The standard constructor, generates an empty object
    RHSSrcPts();

    //! The copy constructor
    RHSSrcPts( const RHSSrcPts &other );

    //! The destructor
    virtual ~RHSSrcPts();

    //! The assignment operator
    RHSSrcPts &operator=( const RHSSrcPts &other );

    //! The copy function to return a pointer to a cop of the object
    virtual RHS *copy() const;

    //! A function that clears out all source points
    void clearSrcPts();
    
    //! A function that sets the accuracy of approximating the point sources
    void setEpsilon( Double epsilon );

    //! A function that sets the vector of point source locations
    void setSrcPts( const vector<DoubleVector> &sourcPts );

    //! A function that adds a point source location
    void addSrcPt( const DoubleVector &src );

    //! A function that returns all currently saved point source locations
    void getSrcPts( vector<DoubleVector> &res ) const;

    //! A function that returns the current accuracy of approximation
    Double getEpsilon() const;

    //! A functiont that evaluates the right hand side at a point x
    Complex eval( DoubleVector x ) const;

  private:
    //! An auxilliary function that approximates a point source by an exponential function
    /*!
     * \p x0 The point source location
     * \p x The evaluation point
     */
    Double deltaFunction (const DoubleVector &x0, const DoubleVector &x ) const;
};

#endif
