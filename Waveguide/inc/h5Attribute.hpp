#ifndef H5_ATTRIBUTE_H
#define H5_ATTRIBUTE_H
#include <H5Cpp.h>

struct H5Attr
{
  H5Attr():name(""), value(0.0), dtype(H5::PredType::NATIVE_DOUBLE){};
  std::string name;
  double value;
  H5::PredType dtype;
};

H5Attr makeAttr( const char* name, double value );
H5Attr maktAttr( const char* name, int value );

#endif
