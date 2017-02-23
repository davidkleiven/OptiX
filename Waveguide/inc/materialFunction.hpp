#ifndef MATERIAL_FUNCTION_H
#define MATERIAL_FUNCTION_H

class MaterialFunction
{
public:
  virtual void getXrayMatProp( double x, double y, double z, double &delta, double &beta ) const = 0;
};
#endif
