#ifndef PML_H
#define PML_H

#include "typedef.h"

class AssemblerHelmholtz;

//! The abstract class for a PML implementation
/*!
 * This class is overloaded by specific PML implementations. The abstract class is use in other parts of the code (e.g. a domain decomposition).
 */
class PML
{
  protected:
    //! The pointer to the parent assembly routine of the 
    const AssemblerHelmholtz *mParent;

    //! The range of the computational domain in x-direction
    vector<DoubleVector> mCompDomain;
    
    //! The range of the interior domain (domain without PML) in x-direction
    vector<DoubleVector> mIntDomain;

  private:
    //! The physical width of the PML
    Double mPMLWidth;

  public:
    //! The standard constructor
    PML();
    
    //! The copy constructor
    PML( const PML &other );

    //! The assignment operator
    PML &operator=( const PML &other );

    //! The destructor
    virtual ~PML() {}

    //! The copy function
    /*!
     * This function copies the (derived) object properly and returns a pointer to
     * the abstract class. This is important to be able to store different types
     * of PML implementations in other parts of the code while not changing the code
     * and consequently makes the code more modular.
     */
    virtual PML *copy() const = 0;

    //! The function to set the parent assembly routine
    void setParent( const AssemblerHelmholtz *parent );

    //! The funciton that sets the physical width of the PML
    virtual void setWidth( Double pmlWidth );

    //! The function that returns the physical width of the PML
    virtual Double getWidth() const;

    // The function that sets the interior domain of the PML
    /*!
     * The input domain should only be the interior domain,
     * the computational domain is computed automatically
     * using mPMLWidth.
     */
    virtual void setDomain( const vector<DoubleVector> &domain );
    
    //! The function that returns the interior domain of the PML
    virtual vector<DoubleVector> getDomain() const;

    //! The function that returns the coefficient for the
    //  diffusion (second derivative) term that needs to be added due to the PML
    //  \param x The evaluation point
    //  \return The vector of values that each direction has to be multiplied with
    virtual ComplexVector getDiffusionCoefficient( DoubleVector x ) const = 0;
    
    //! The function that returns the coefficient for the
    //  convection (first derivative) term that needs to be added due to the PML
    //  \param x The evaluation point
    //  \return The vector of values that each direction has to be multiplied with
    virtual ComplexVector getConvectionCoefficient( DoubleVector x ) const = 0;
    
    //! The function that returns the coefficient for the
    //  reaction (no derivative) term that needs to be added due to the PML
    //  \param x The evaluation point
    //  \return The the value that the reaction term has to be multiplied with
    virtual Complex getReactionCoefficient( DoubleVector x ) const = 0;
    
    //! The function that returns the coefficient for the
    //  source function (rhs) term that needs to be added due to the PML
    //  \param x The evaluation point
    //  \return The the value that the source term has to be multiplied with
    virtual Complex getSourceCoefficient( DoubleVector x ) const = 0;
};

#endif

