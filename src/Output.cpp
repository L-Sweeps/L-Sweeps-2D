#include "Output.h"

// The standard constructor, generates an empty object
Output::Output(): mDD(NULL)
{

}

// The copy constructor, copies all memeber variables
Output::Output( const Output &other ): mDD(other.mDD), mFilename(other.mFilename)
{

}

// The destructor, nothing to be done
Output::~Output()
{

}

// The assignment operator, copies all member variables
Output &Output::operator=( const Output &other )
{
  // only copy values if other is not equal to this object
  if( this != &other )
  {
    // copy the member variables
    mDD = other.mDD;
    mFilename = other.mFilename;
  }
  return *this;
}

// sets the DD that is used for this object
void Output::setDD( const DD *DD )
{
  mDD = DD;
}

// sets the filename for the output file
void Output::setFilename( string name )
{
    mFilename = name;
}

// returns the currently set filename
string Output::getFilename() const
{
    return mFilename;
}
