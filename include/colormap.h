#ifndef COLORMAP_H
#define COLORMAP_H

#include <glm/glm.hpp>
using namespace glm;

class Colormap{
private:
  int nsamples;
  double step;
  float *domain;
  float gamma;
  vec4  *range;
  vec4 computecolor(float x);
public:
  Colormap(float gamma);
  vec4 colorof(double x);
  void destroy();
};

#endif