#include <math.h>
#include "deconvolve.h"
#include <cstdio>
#include <glm/glm.hpp>
#define pi (3.141592653589793)
#define gaus(x,stdv) (exp(-((x)*(x))/(2*(stdv)*(stdv)))/((stdv)*sqrt(2*pi)))

Imdim::Imdim(int a0, int a1, int a2) : a0(a0), a1(a1), a2(a2), w0(1), w1(a0), w2(a0*a1), w3(a0*a1*a2){
  // #
}
float PSF::eval(float x, float y, float z, float dx, float dy, float dz){
  // simply gaussian deconvolution with sigma=2.
  // float d = sqrt(dx*dx + dy*dy + dz*dz);
  // return gaus(d, 2);
  glm::vec3 p(dx,dy,dz);
  return exp(-0.5 * glm::dot(p,(shape*p)));
}

PSF::PSF() : shape(){

}

PSF::PSF(glm::mat3x3 s) : shape(s){

}

inline int max(int x, int y){
  return (x<y)?y:x;
}
inline int min(int x, int y){
  return (x<y)?x:y;
}
float PSF::convolve(float *data, float x, float y, float z, Imdim dim, bool verb){
  int xx = x;
  int yy = y;
  int zz = z;
  // printf("%d %d %d\n", dim.w0, dim.w1, dim.w2);
  // return data[xx*dim.w0 + yy*dim.w1 + zz*dim.w2];

  int R = 3;

  int xmin=max(x-R,0);
  int xmax=min(x+R,dim.a0-1);

  int ymin=max(y-R,0);
  int ymax=min(y+R,dim.a1-1);

  int zmin=max(z-R,0);
  int zmax=min(z+R,dim.a2-1);

  float val = 0;
  float sum = 0;
  // printf("%d %d; %d %d; %d %d\n", xmin, xmax, ymin, ymax, zmin, zmax);
  for(int xx=xmin; xx<=xmax; xx++){
    for(int yy=ymin; yy<=ymax; yy++){
      for(int zz=zmin; zz<=zmax; zz++){
        float k = eval(x,y,z,x-xx,y-yy,z-zz);
        val += k * data[xx*dim.w0 + yy*dim.w1 + zz*dim.w2];
        // if(verb)printf("  val += %.5f * %.5f\n", k, data[xx*dim.w0 + yy*dim.w1 + zz*dim.w2]);
        sum += k;
        // val = data[dim.w0*(int(x) + xx) + dim.w1*(int(y) + yy) + dim.w2*(int(z) + zz)];
        // sum = 1;
        // printf("%d / \n", dim.w0*(int(x) + xx) + dim.w1*(int(y) + yy) + dim.w2*(int(z) + zz));
        // sum += 1;
        // float k = eval(x,y,z,x-xx,y-yy,z-zz);
        // sum += k;
        // val += k * data[dim.w0*(int(x) + xx) + dim.w1*(int(y) + yy) + dim.w2*(int(z) + zz)];
      }
    }
  }
  if(!sum)sum=1;
  // if(verb)printf("%.5f / %.5f = %.5f\n", val, sum, val/sum);
  val /= sum;
  return val;
}

#undef pi
#undef gaus