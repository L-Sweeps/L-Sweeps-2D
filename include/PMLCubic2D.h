#ifndef PMLQUADRATIC_H
#define PMLQUADRATIC_H

#include<vector>
#include<map>
#include "PML.h"
#include "typedef.h"
#include "AssemblerHelmholtz.h"

using namespace std;

//! A class that implements the PML in 2D using a quadratic profile
class PMLCubic2D: public PML
{
  protected:
    //! The absorption constant
    Double mC;

  public:
    //! The standard constrctor, generates an empty object
    PMLCubic2D() {}

    //! The copy constructor
    PMLCubic2D( const PMLCubic2D &other );

    //! The destructor
    virtual ~PMLCubic2D() {}

    //! The copy function
    /*!
     * This function generates a copy of this object and returns a pointer to it.
     */
    virtual PML *copy() const;
    
    //! The assignment operator
    PMLCubic2D &operator=( const PMLCubic2D &other );

    //! A function that sets the absorption coefficient
    /*!
     * \param C The absorption constant to be set
     */
    void setAbsorptionCoeff( Double C );
  
    //! A function that returns the absorption coefficient
    /*!
     * \return The absorption coefficient
     */
    Double getAbsorptionCoeff() const { return mC; }
 
    //! The function that returns the diffusion coefficient
    /*!
     * \sa PML
     */
    virtual ComplexVector getDiffusionCoefficient( DoubleVector x ) const;

    //! The function that returns the convection coefficient
    /*!
     * \sa PML
     */
    virtual ComplexVector getConvectionCoefficient( DoubleVector x ) const;

    //! The function that returns the reaction coefficient
    /*!
     * \sa PML
     */
    virtual Complex getReactionCoefficient( DoubleVector x ) const;

    //! The function that returns the source coefficient
    /*!
     * \sa PML
     */
    virtual Complex getSourceCoefficient( DoubleVector x ) const;

  private:
    //! A function that evaluates \f$\alpha_1\f$ according to the notation
    //  used [here](https://dspace.mit.edu/openaccess-disseminate/1721.1/115980)
    /*!
     * \param x The evaluation point
     * \return The value of \f$\alpha_1\f$
     */
    Complex alpha1( DoubleVector x ) const;

    //! A function that evaluates the derivative of \f$\alpha_1\f$ according to 
    //  the notation used [here](https://dspace.mit.edu/openaccess-disseminate/1721.1/115980)
    /*!
     * \param x The evaluation point
     * \return The value of the derivative of \f$\alpha_1\f$
     */
    Complex alpha1Dx( DoubleVector x ) const;

    //! A function that evaluates \f$\alpha_2\f$ according to the notation
    //  used [here](https://dspace.mit.edu/openaccess-disseminate/1721.1/115980)
    /*!
     * \param x The evaluation point
     * \return The value of \f$\alpha_2\f$
     */
    Complex alpha2( DoubleVector x ) const;

    //! A function that evaluates the derivative of \f$\alpha_2\f$ according to 
    //  the notation used [here](https://dspace.mit.edu/openaccess-disseminate/1721.1/115980)
    /*!
     * \param x The evaluation point
     * \return The value of the derivative of \f$\alpha_2\f$
     */
    Complex alpha2Dy( DoubleVector x ) const;

    //! A function that evaluates \f$\sigma_1\f$ according to the notation
    //  used [here](https://dspace.mit.edu/openaccess-disseminate/1721.1/115980)
    /*!
     * \param x The evaluation point
     * \return The value of \f$\sigma_1\f$
     */
    Double sigma1( DoubleVector x, Double omega ) const;
    
    //! A function that evaluates the derivative of \f$\sigma_1\f$ according to
    //  the notation used [here](https://dspace.mit.edu/openaccess-disseminate/1721.1/115980)
    /*!
     * \param x The evaluation point
     * \return The value of the derivative of \f$\sigma_1\f$
     */
    Double sigma1Dx( DoubleVector x, Double omega ) const;
    
    //! A function that evaluates \f$\sigma_2\f$ according to the notation
    //  used [here](https://dspace.mit.edu/openaccess-disseminate/1721.1/115980)
    /*!
     * \param x The evaluation point
     * \return The value of \f$\sigma_2\f$
     */
    Double sigma2( DoubleVector x, Double omega ) const;

 
    //! A function that evaluates  the derivative of \f$\sigma_2\f$ according to
    //  the notation used [here](https://dspace.mit.edu/openaccess-disseminate/1721.1/115980)
    /*!
     * \param x The evaluation point
     * \return The value of the derivative of \f$\sigma_2\f$
     */
    Double sigma2Dy( DoubleVector x, Double omega ) const;
};

#endif

