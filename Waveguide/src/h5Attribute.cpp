#include "h5Attribute.hpp"

H5Attr makeAttr( const char* name, double value )
{
  H5Attr attr;
  attr.value = value;
  attr.dtype = H5::PredType::NATIVE_DOUBLE;
  attr.name = name;
  return attr;

  return attr;
}

H5Attr makeAttr( const char* name, int value )
{
  H5Attr attr;
  attr.value = value;
  attr.dtype = H5::PredType::NATIVE_INT;
  attr.name = name;
  return attr;
}
