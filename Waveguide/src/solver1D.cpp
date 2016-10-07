#include "solver1D.hpp"

Solver1D::~Solver1D()
{
  delete solution;
}

void Solver1D::fillJsonObj( Json::Value &obj ) const
{
  obj["name"] = name;
  obj["eigenvalue"] = eigenvalue;
}
