#ifndef FIELD_MONITORS_H
#define FIELD_MONITORS_H
#include <armadillo>
#include <meep.hpp>
#include <string>

class FieldMonitor
{
public:
  FieldMonitor( unsigned int N1, unsigned int N2 );
  ~FieldMonitor();

  /** Returns the intensity */
  arma::mat& get() { return *poynting; };

  /** Set intensity */
  void setIntensity( const meep::fields &field );

  /** Set the displacement vectors in each direction */
  void setDisplacementVectors( const meep::vec &v1, const meep::vec &v2 );

  /** Sets a point in the monitor plane */
  void setOrigin( const meep::vec &orig ){ origin = orig; };

  /** Get name of monitor */
  const std::string& getName() const { return name; };

  /** Set the name of the monitor */
  void setName( const char* newname ){ name = newname; };
private:
  meep::vec displacement1;
  meep::vec displacement2;
  meep::vec origin;
  arma::mat *poynting{NULL};
  std::string name{""};
};
#endif
