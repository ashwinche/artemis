#include <teem/meet.h>
#include "experiment.h"
class Flow{
  ArExperiment *exp;
  Nrrd *buff;
  Nrrd *compute_flow(int low, int high);
};