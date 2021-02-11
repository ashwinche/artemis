#ifndef DECONV
#define DECONV

#include <glm/glm.hpp>
struct Imdim{
  Imdim(int a0, int a1, int a2);
  int a0,a1,a2,w0,w1,w2,w3;
};
class Deconvolve{
  void deconvolve(float *data);
};

struct PSF{
  PSF();
  PSF(glm::mat3x3);
  glm::mat3 shape;
  float eval(float x, float y, float z, float dx, float dy, float dz);
  float convolve(float *data, float x, float y, float z, Imdim dim, bool verb=false);
};

#endif