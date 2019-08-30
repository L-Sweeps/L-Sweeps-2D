#ifndef ASSEMBLERHELMHOLTZ_H
#define ASSEMBLERHELMHOLTZ_H

#include "PML.h"
#include "typedef.h"
#include "Assembler.h"
#include "Function.h"

//! An abstract class for the assembly of a Helmholtz problem
/*! 
 * This class assumes that the problem is solved with PML absorbing boundary
 * conditions around the entire domain and assembles the PDE
 * \f{equation}{ -\Delta u -m\omega^2u=f\quad\text{in }\Omega.\f}
 * where \f$\Omega\f$ is a square.
 * Note that this class is still abstract, as it does not define all pure
 * virtual functions in Assembler.
 */
class AssemblerHelmholtz : public Assembler
{
  protected:
    //! The object realizing the PML absorbing boundary conditions
    PML *mPML;

    //! The function \f$m\f$ as a FunctionPtr
    Function* mM;

    //! The frequency \f$\omega\f$
    Double mOmega;

  public:
    //! The standard constructor
    /*! 
     * This constructor generates an empty object. Properties and values
     * have to be set with the set functions.
     */
    AssemblerHelmholtz();

    //! The copy constructor
    /*!
     * Copies all member variables correctly.
     */
    AssemblerHelmholtz( const AssemblerHelmholtz &other );

    //! The destructor
    /*!
     * Properly frees up of all allocated memery before destruction.
     */
    virtual ~AssemblerHelmholtz();

    //! The assignment operator
    /*!
     * Properly copies \p other and assigns it to this object.
     */
    AssemblerHelmholtz &operator=( const AssemblerHelmholtz &other );
  
    //! The function to set the PML in this object
    /*!
     * \param pml The PML to be set. This function properly copies the dereferenced pointer.
     */
    void setPML( PML *pml );

    //! The function to set the function \f$m\f$
    /*!
     * \param m The pointer to the function.
     */
    void setM( Function* m ); 
    
    //! The function to set the frequency \f$\omega\f$
    /*!
     * \param omega The value for the frequency
     */
    void setOmega( Double omega );
 
    //! The function to return the function \f$m\f$
    /*!
     * \return the FunctionPtr for \f$m\f$
     */
    Function* getM() const { return mM; }
    
    //! The function to return the frequency \f$\omega\f$
    Double getOmega() const { return mOmega; }
};

#endif
