#ifndef FAR_FIELD_PARAMETERS_H
#define FAR_FIELD_PARAMETERS_H

/** Parameters used when storing the far field */
struct FarFieldParameters
{
  double padValue{0.0};
  double phiMin{-180.0};
  double phiMax{180.0};
};
#endif
