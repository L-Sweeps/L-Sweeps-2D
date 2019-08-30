#ifndef OUTPUT_H
#define OUTPUT_H

#include "DD.h"
#include <string>

using namespace::std;


//! The abstract class for an output routine
/*!
 * This class is overloaded by specific output routines.
 */
class Output
{
  protected:
    //! The pointer to the parent domain decomposition
    const DD *mDD;
    
    //! The filename used for output
    string mFilename;

  public:
    //! The standard constructor
    Output();

    //! The copy constructor
    Output( const Output &other );

    //! The destructor
    ~Output();

    //! The assignment operator
    Output &operator=( const Output &other );

    //! The copy function
    /*!
     * This function is written so that objects of derived
     * classes are properly copied.
     */
    virtual Output *copy() const = 0;

    //! The function that sets the parent domain decomposition
    void setDD( const DD *DD );

    //! The function that sets the filename
    void setFilename( string name );

    //! The functino that returns the filename
    string getFilename() const;

    //! The function that prints the output to file
    virtual void print() const = 0;

    //! The function that assigns values to the output
    /*!
     * The set values are used to generate the output.
     */
    virtual void setVals( const vector<ComplexVector> &vals ) const = 0;
};

#endif
