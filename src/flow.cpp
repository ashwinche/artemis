#include "flow.h"

Nrrd* Flow::compute_flow(int low, int high){
  if(low<0)low = 0;
  if(high>exp->high)high = exp->high;

  buff = exp->copy(0);
  return 0;
  // Filter filter;
}