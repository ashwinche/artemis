#include "view.h"
#include "colormap.h"

Colormap::Colormap(float gamma) : gamma(gamma){
  nsamples = 512;
  domain = new float[nsamples];
  range  = new vec4[nsamples];
  step = 2.0/double(nsamples);
  for(int i=0;i<nsamples;++i){
    float x = (float(i)/float(nsamples));
    range[i] = computecolor(x*2.0);
  }
}
vec4 Colormap::colorof(double x){
  if(x>=6.f){
    x -= 6.f;
    int n = int(x/step);
    if(n>=nsamples)n=nsamples-1;
    return vec4(range[n].x, range[n].z, range[n].y, range[n].w);
  }
  else if(x>=4.f){
    x -= 4.f;
    int n = int(x/step);
    if(n>=nsamples)n=nsamples-1;
    return vec4(range[n].z, range[n].y, range[n].x, range[n].w);
  }
  else if(x>=2.f){
    x -= 2.f;
    int n = int(x/step);
    if(n>=nsamples)n=nsamples-1;
    return vec4(range[n].y, range[n].x, range[n].z, range[n].w);
  }
  else{
    int n = int(x/step);
    if(n<0)n=0;
    if(n>=nsamples)n=nsamples-1;
    return range[n];
  }
}
static float sq(float x){
  return x*x;
}
vec4 Colormap::computecolor(float x){
  vec4 color(0,0,0,0);
  if(x<=1.f){
    color.x = x*x;
    color.y = x;
    color.z = 1.f - sqrt(x);
    float w = sqrt(x);
    color.w = pow(w,gamma);
  }
  else{
    x -= 1.f;
    color.x = sqrt(fmax(1-fabs(2*x),0));
    color.y = sqrt(fmax(1-fabs(2*(x-0.5f)),0));
    color.z = sqrt(fmax(1-fabs(2*(x-1.f)),0));
    color.w = 0.7f;
    // color=vec4(1,1,1,1);
  }
  return color;
}
void Colormap::destroy(){
  delete[] domain;
  delete[] range;
  domain=0;
  range=0;
}