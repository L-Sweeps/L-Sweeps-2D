#ifndef OUTPUTVTK_H
#define OUTPUTVTK_H

#include "Output.h"

//! A class that generates VTK output for Paraview
class OutputVTK: public Output
{
  protected:
    //! The (solution) values used for the output
    /*!
     * This vector stores the values corresponding to the 
     * degrees of freedom in each subdomain. The outer vector
     * is a vector over all subdomains.
     */
    mutable vector<ComplexVector> mVals;

  public:
    //! The standard constructor
    /*!
     * Constructs an empty object.
     */
    OutputVTK();

    //! The copy constructor
    OutputVTK( const OutputVTK &other );

    //! The destructor
    ~OutputVTK();

    //! The assignment operator
    OutputVTK &operator=( const OutputVTK &other );

    //! This function properly copies an object and returns a pointer to it
    Output *copy() const;

    //! The function that prints the output file
    /*!
     * This function uses the filename and assigned values
     * to generate an VTK file to be used in Paraview and
     * allows for the visualization of a solution.
     */
    void print() const;

    //! The function that sets the values
    /*!
     * \p vals The values that should be set. The outer vector is over all subdomains.
     */
    void setVals( const vector<ComplexVector> &vals ) const;
    
    //! The function that returns the currently saved values
    /*!
     * \return The currently stored values in mVals
     */
    vector<ComplexVector> getVals() const;

  private:
    //! This is an auxilliary function to generate files for 2D DDs
    void print2D() const;
};

#endif
