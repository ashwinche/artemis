#include "blob.h"
#include "filter.h"
#include <teem/meet.h>


class Estimator{
public:
  ScaleBlob fit(Nrrd* source, ScaleBlob* in);
  void select_scales(Nrrd* source, ScaleBlob* in, ArFilter &filter);
};