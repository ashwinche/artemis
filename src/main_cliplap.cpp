#include <cstdio>
#include <cstdlib>
#include <teem/meet.h>
#include <cmath>
#include "filter.h"

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
    maxv = fmax(maxv, out[j]);
  }

  // normalize max value to 1.0
  for(int j=0;j<a1*a2*a3;++j){
    out[j] /= maxv;
  }  

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

  float maxv=0;
  for(int i=0;i<a0*a1*a2;++i){
    out[i] = (float)(in[i] / 32768.0);
    maxv = fmax(maxv, out[i]);
  }

  // normalize max value to 1.0
  for(int j=0;j<a0*a1*a2;++j){
    out[j] /= maxv;
  }  

  nrrdWrap_va(nout, out, nrrdTypeFloat, 3, a0, a1, a2);
  nrrdNuke(nin);
  return nout;
}

static Nrrd* load_nrrd(const char *filename){
  // printf("load %s\n", filename);
  Nrrd *nrrd = nrrdNew();
  if (nrrdLoad(nrrd, filename, NULL)) {
    char *err = biffGetDone(NRRD);
    fprintf(stderr, "error reading \"%s\":\n%s", filename, err);
    free(err);
    exit(0);
  }
  int dim = nrrd->dim;
  if(dim == 4){
    return convert_short4_to_float3_and_destroy(nrrd);
  }else if(dim == 3){
    return convert_short3_to_float3_and_destroy(nrrd);
  }else{
    printf("unable to read nrrd of dim=%d. exit.\n", dim);
    exit(0);
  }
}
static void emit_nrrd(Nrrd *nrrd, char *file){
  Nrrd *nrrdfloat = nrrd;
  Nrrd *nrrdshort = nrrdNew();
  nrrdConvert(nrrdshort, nrrdfloat, nrrdTypeUShort);
  int size = nrrdshort->axis[0].size * nrrdshort->axis[1].size * nrrdshort->axis[2].size;
  short *datashort = (short*)nrrdshort->data;
  float *datafloat = (float*)nrrdfloat->data;
  for(int i=0;i<size;++i){
    datashort[i] = (unsigned short) (datafloat[i] * 32767.0);
  }
  // nrrdCopy(out, asfloat);
  nrrdSave(file, nrrdshort, NULL);
  nrrdNuke(nrrdshort);
}
int main(int argc, char** argv){
  setvbuf(stdout, NULL, _IONBF, 0);  
  // printf("hello.");
  // exit(0);
  if(argc != 4){
    printf("USAGE: ./cliplap [INPUT] [OUTPUT] [SCALE]\n");
    exit(0);
  }
  char *filein  = argv[1];
  char *fileout = argv[2];

  float scale = atof(argv[3]);
  // printf("load_nrrd..");
  Nrrd *nin = load_nrrd(filein);
  // printf("nrrdNew..");

  ArFilter filter;
  Nrrd *base = nrrdNew();
  // printf("nrrdCopy..");
  nrrdCopy(base, nin);
  // printf("filter.init..");
  filter.init(base);

  // printf("done..");

  if(scale > 0){
    filter.max1();
    DiscreteKernel kernel = filter.gaussian(scale, int(scale*4));
    filter.set_kernel(kernel);
    filter.filter();
    kernel.destroy();
    filter.laplacian3d();
    printf("%.4f", scale);
  }
  else{
    int last_np  = 0;
    int last_dnp = 0;
    int last_ddnp = 0;
    for(float f = 0.5f; f< 15.f; f*= 1.33f){
      // printf("capture..");
      filter.capture(nin);
      DiscreteKernel kernel = filter.gaussian(f, int(f*4));

      // printf("median..");
      filter.max1();    
      filter.set_kernel(kernel);
      // printf("filter..");
      filter.filter();
      // printf("laplacian..");
      kernel.destroy();
      filter.laplacian3d();
      // printf("count..");
      int np = filter.count_blobs();

      int dnp = last_np - np; 
      int ddnp = last_dnp - dnp;
      int dddnp = last_ddnp - ddnp;

      // printf("%.2f\t%d\t%d\t%d\t%d\n", f, np, dnp, ddnp, dddnp);
      printf("%.2f\t%d\n", f, np);
      if(dnp > 0){
        // printf("%.4f", f);
        // printf(" **** \n");
        // break;
      }
      last_np = np;
      last_dnp = dnp;
      last_ddnp = ddnp;
    }
  }

  filter.normalize();
  filter.commit(nin);

  emit_nrrd(nin, fileout);
}