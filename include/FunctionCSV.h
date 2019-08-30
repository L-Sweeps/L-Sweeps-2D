#ifndef FUNCTIONCSV_H
#define FUNCTIONCSV_H

#include "Function.h"

//! A class that realizes a function defined by a CSV file
/*!
 * The CSV file is read in as a matrix and a uniform mesh of squares.
 * The domain of the mesh is given with the file when the matrix is
 * generated. A piecewise constant function is defined on this mesh
 * where the values are inherited from the matrix.
 *
 * In addition, a domain of definition is defined. When the function
 * is evaluated for a value in the domain of definition, it depends
 * if this value is inside or outside the defined mesh. If it is inside,
 * the corresponding value is returned. If it is outside, the closest point
 * in the mesh is found in the direction normal to the mesh. Then the
 * corresponding value to this point is returned. For a value outside the
 * domain of defintion, the closest point in the domain of definition
 * in the normal direction is found and the corresponding value returned.
 */
class FunctionCSV : public Function
{
  protected:
    //! The filename used for the CSV file
    string mFilename;

    //! The domain read in with the file and used to define the mesh
    vector<DoubleVector> mFileDom;

    //! The number of values in the file
    UIntVector mFileN;

    //! The mesh size in each direction of the mesh
    DoubleVector mH;

  private:
    //! The domain of definition
    mutable vector<DoubleVector> mDom;

    //! The subset of the mesh covering the domain of definition
    mutable vector<DoubleVector> mDomA;

    //! The matrix holding the values in the mesh
    mutable ComplexMatrix mA;

    //! The minimum real value in the matrix
    mutable Double mMinReal;
    
    //! The maximum real value in the matrix
    mutable Double mMaxReal;
    
    //! The minimum imaginary value in the matrix
    mutable Double mMinImag;
    
    //! The maximum imaginary value in the matrix
    mutable Double mMaxImag;
    
    //! The flag determining if the domain of definition is set 
    mutable bool mDomSet;

  public:
    //! The standard constructor
    /*!
     * This function sets up an empty object.
     */
    FunctionCSV();

    //! The copy constructor
    FunctionCSV( const FunctionCSV &other );
    
    //! The destructor
    virtual ~FunctionCSV();

    //! The assignment operator
    FunctionCSV &operator=( const FunctionCSV &other );

    //! The function that clears the cache
    virtual void clearCache() const;

    //! The function that sets the file
    /*!
     * This function prepares the file to be read in and
     * sets up the mesh in the given domain.
     *
     * \p filename The filename including the .csv extension
     * \p fileDom The domain associated with the file
     */
    void setFile( string filename, const vector<DoubleVector> &fileDom );

    //! The function return the minimum real value
    virtual Double getMinReal() const;
    
    //! The function return the maximum real value
    virtual Double getMaxReal() const;
    
    //! The function return the minimum imaginary value
    virtual Double getMinImag() const;
    
    //! The function return the maximum imaginary value
    virtual Double getMaxImag() const;

    //! The function that sets the domain of definition
    /*!
     * This function reads the required part of the matrix from file
     * and sets up all required information.
     *
     * \p dom The domain to be set
     */
    virtual void setDomain( const vector<DoubleVector> &dom ) const;
    
    //! The function returning the domain of definition
    /*!
     * \p res The result
     */
    virtual void getDomain( vector<DoubleVector> &res ) const;

    //! The copy function to properly copy this object
    virtual Function *copy() const;

    //! The evaluation function
    /*!
     * \p x The point at which the function should be evaluated
     */
    virtual Complex eval( DoubleVector x ) const;
};

#endif
