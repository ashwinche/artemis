#include "experiment.h"
#include <SFML/Graphics.hpp>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>

static Nrrd* convert_short4_to_float3_and_destroy(Nrrd *nin){
  // printf("hi!\n");
  int a0 = nin->axis[0].size;
  int a1 = nin->axis[1].size;
  int a2 = nin->axis[2].size;
  int a3 = nin->axis[3].size;

  Nrrd *nout = nrrdNew();

  short *in  = (short*)nin->data;
  float *out = new float[a1*a2*a3];

  float maxv=0;
  for(int j=0;j<a1*a2*a3;++j){
    out[j] = (float)(in[j*a0+0] / 32768.0);
    maxv = max(maxv, out[j]);
  }

  // normalize max value to 1.0
  for(int j=0;j<a1*a2*a3;++j){
    out[j] /= maxv;
  }  

  // printf("AWENRKAWEJRKAWEKR");
  // for(int z=0;z<a2;z++){
  //   float sum = 0;
  //   for(int x=0;x<a0;x++){
  //     for(int y=0;y<a1;y++){
  //       // sum += out[x + y*a0 + z*a0*a1];
  //       sum = max(sum, out[x + y*a0 + z*a0*a1]);
  //     }
  //   }
  //   for(int x=0;x<a0;x++){
  //     for(int y=0;y<a1;y++){
  //       out[x + y*a0 + z*a0*a1] /= sum;
  //     }
  //   }
  //   printf("%.2f\n", sum);
  // }

  nrrdWrap_va(nout, out, nrrdTypeFloat, 3, a1, a2, a3);
  // printf("convert %d %d %d\n", nout->axis[0].size, nout->axis[1].size, nout->axis[2].size);
  nrrdNuke(nin);
  // printf("created %p.\n", nout);
  return nout;
}

static Nrrd* convert_short3_to_float3_and_destroy(Nrrd *nin){
  // printf("hello3!\n");
  int a0 = nin->axis[0].size;
  int a1 = nin->axis[1].size;
  int a2 = nin->axis[2].size;

  Nrrd *nout = nrrdNew();

  short *in  = (short*)nin->data;
  float *out = new float[a0*a1*a2];
  // printf("yellow");
  float maxv=0;
  for(int i=0;i<a0*a1*a2;++i){
    out[i] = (float)(in[i] / 32768.0);
    maxv = max(maxv, out[i]);
  }
  // printf("hi there");

  // normalize max value to 1.0
  for(int j=0;j<a0*a1*a2;++j){
    out[j] /= maxv;
  }

  // printf("again.");
  nrrdWrap_va(nout, out, nrrdTypeFloat, 3, a0, a1, a2);
  nrrdNuke(nin);
  return nout;
}

static Nrrd* load_nrrd(const char *filename){
  printf("load %s\n", filename);
  Nrrd *nrrd = nrrdNew();
  if (nrrdLoad(nrrd, filename, NULL)) {
    char *err = biffGetDone(NRRD);
    fprintf(stderr, "error reading \"%s\":\n%s", filename, err);
    free(err);
    exit(0);
  }
  int dim = nrrd->dim;
  // printf("hello");
  if(dim == 4){
    // printf("5");
    return convert_short4_to_float3_and_destroy(nrrd);
  }else if(dim == 3){
    // printf("4");
    return convert_short3_to_float3_and_destroy(nrrd);
  }else{
    printf("unable to handle nrrd of dim=%d. exit.\n", dim);
    exit(0);
  }
  // printf("loaded");
}


/** Returns length of first sequence of 'c' in str.
  * Eg. f("aabbcccdefccccc", 'c') = 3.
  */ 
static int length_first_repeated_sequence(const char *str, const char c){
  int digits = 0;
  const char *p = str;
  while(*p){          // while not hit null-terminator.
    if(*p == c){      // count number of question marks.    
      ++digits;
    }
    else if(digits){  // no longer a question mark, after
      return digits;  // we started counting (because digits > 0).
    }                 // return digits.
    ++p;
  }
  return digits;
}

std::string resolve_filename(std::string path, int index){
  int digits = length_first_repeated_sequence(path.c_str(),'?');   // number of digits;
char counter[digits+1];
  memset(counter,'\0',digits);      // \0 \0 .... digits .... \0
  snprintf(counter, sizeof(counter), "%0*d", digits, index);
  std::string filename = path;
  int p;
  while((p=filename.find('?')) != std::string::npos){
    filename = filename.replace(p,digits,counter);
  }
  return filename;
}

ArExperiment::ArExperiment(std::string path, int low, int high, int mem_cap){
  paths     = new std::string[high-low+1];
  for(int i=low;i<=high;++i){
    paths[i-low] = resolve_filename(path, i);
    // printf("path = %s\n", paths[i-low].c_str());
  }

  frames = new NrrdFrame[mem_cap];
  for(int i=0;i<mem_cap;++i){
    frames[i].n        = -1;
    frames[i].accessed = 0;
    frames[i].path     = "";
    frames[i].nrrd     = 0;
  }
  this->low     = low;
  this->high     = high;
  this->nframes = mem_cap;
  this->npaths  = (high-low+1);
  this->time    = 0;
  this->filepath = path;
}

void ArExperiment::interpolate(Nrrd* f0, Nrrd *f1, Nrrd *fx, float alpha){
  using namespace glm;
  float *d0 = (float*)f0->data;
  float *d1 = (float*)f1->data;
  float *dx = (float*)fx->data;

  int a0 = f0->axis[0].size;
  int a1 = f0->axis[1].size;
  int a2 = f0->axis[2].size;

  int w2 = a1*a0;
  int w1 = a0;

  for(int i=0;i<a0*a1*a2;i++)dx[i] = 0;
  for(int i=0;i<a0*a1*a2;i++)dx[i] = alpha*d1[i] + (1-alpha)*d0[i];
  
  // for(int x=0;x<a0;x++){
  //   for(int y=0;y<a1;y++){
  //     for(int z=0;z<a2;z++){
  //       vec3 v(10,0,0);
  //       int xx = int(x + v.x);
  //       int yy = int(y + v.y);
  //       int zz = int(z + v.z);
  //       if(xx>=0 && yy>=0 && zz>=0 && xx<a0 && yy<a1 && zz<a2){
  //         fx[zz*w2 + yy*w1 + zz] += 
  //       }
  //     }
  //   }
  //   // dx[i] = alpha*d1[i] + (1-alpha)*d0[i];
  // }
}

Nrrd* ArExperiment::get(float n, bool force){

  if(n - int(n)   < 1e-3)n=int(n);
  if(int(n+1) - n < 1e-3)n = int(n+1);
  // printf("get frame %d\n",n);
  ++time;
  int min_i   = 0;
  int min_acc = frames[0].accessed;
  // look for frame n in memory.
  for(int i=0;i<nframes;++i){
    if(frames[i].accessed < min_acc){
      min_i = i;
      min_acc = frames[i].accessed;
    }
    if(frames[i].n == n){
      frames[i].accessed = time;
      if(force){    // found it, but re-load it anyway.
        min_i = i;
        break;
      }
      return frames[i].nrrd;
    }
  }
  // frame n not in memory. load from disk.
  int i=min_i;
  if(frames[i].nrrd){
    nrrdNuke(frames[i].nrrd);
  }
  // printf("loading %s\n",paths[n-low]);
  frames[i].n = n;
  frames[i].accessed = time;
  frames[i].path = paths[int(n)-low];
  frames[i].nrrd = load_nrrd(paths[int(n)-low].c_str());

  if(n != int(n)){
    // fractional time. perform interpolation.
    float alpha = n-int(n);

    Nrrd* f0 = load_nrrd(paths[int(n)-low  ].c_str());
    Nrrd* f1 = load_nrrd(paths[int(n)-low+1].c_str());
    Nrrd* fx = frames[i].nrrd;

    interpolate(f0,f1,fx, alpha);

    nrrdNuke(f0);
    nrrdNuke(f1);
  }

  filter.init(frames[i].nrrd);
  // filter.threshold(0.05,1.0);
  // filter.threshold(0.1,1.0);
  filter.commit();

  // printf("loaded frame %d: %u, %s, %p\n", n, time, paths[n-low].c_str(), frames[i].nrrd);
  return frames[i].nrrd;
}

std::string ArExperiment::getfilepath(int n){
  // printf("path = %s\n", paths[0].c_str());
  return paths[n-low];
}
// std::string ArExperiment::gettgmmpath(int n){
//   return tgmmpaths[n-low];
// }
// std::string ArExperiment::getfilenumber(int n){
//   return numbers[n-low];
// }

Nrrd* ArExperiment::copy(int n){
  Nrrd *src = get(n);
  Nrrd *dst = nrrdNew();
  nrrdCopy(dst, src);

  return dst;
}