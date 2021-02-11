#include <cmath>
#include <teem/nrrd.h>
#include <glm/glm.hpp>
#include <vector>
#include <unordered_set>
#include "blob.h"
#include "bsptree.h"
#include "deconvolve.h"


#ifndef FILTER_H
#define FILTER_H

class DiscreteKernel{
public:
  DiscreteKernel();
  ~DiscreteKernel();
  int radius;
  int support;
  double *data;
  double *temp;
  void destroy();
};


class ArFilter{
public:
  struct{
    Nrrd  **nrrd;
    float **buff;   // transient internal buffers for computing
    int     nbuf;   // successive image filters. There are no
    int     curr;   // guarantees about the contents of these images,
                    // except that buff[curr] holds the result of the
                    // most recent filter.
    
    DiscreteKernel kernel;

    NrrdAxisInfo *a;
    int w0,w1,w2,w3,w4;
    int a0,a1,a2,a3;
    bool alive;
  }self;
  int itempbuf(int c);
  int itempbuf();
  double comp_max_laplacian(float *data);

  void label_blobs(float *data, int* labelled);
  void label_blobs_lindeberg(float *data, int* labelled);
public:
  ArFilter();
  void conv2d(float *nin, float *nout, int xmax, int ymax, int zmax, int xstep, int ystep, int zstep, DiscreteKernel kernel);
  void convolve(Nrrd *nin, Nrrd *nout, DiscreteKernel kernel);
  DiscreteKernel gaussian(double sigma, int radius, int d=0);
  DiscreteKernel interpolation();
  void set_kernel(DiscreteKernel k);
  void filter();

  void posterize(int nlevels);
  void normalize(double power=1.0);
  void scale(float s);
  void threshold(float min, float max);
  void positive(int channel=0);
  void negative(int channel=0);
  void binary(int channel=0);
  void laplacian3d(int boundary = 0);
  void pass();
  void hessian3d(int boundary = 1);
  void laplacianmasked(float scale = 0);
  void max1();
  void median1();
  void maxima();
  void print();
  void clear();
  void show_blobs(int mode=0);

  void deconvolve();


  void difference_image(Nrrd* x);

  int count_connected_components();
  std::vector<std::unordered_set<int>> label_connected_levelsets(float *data, int *labels);
  int count_blobs();


  // void gaussian(float scale);
  // void lapofgaussian(float scale);
  void lapofgaussian_masked(float sigma, int multiply);

  // raster operations
  void rasterlineadd(vec3 a, vec3 b, float va, float vb);

  std::vector<glm::ivec3> find_nonzero();
  std::vector<glm::ivec3> find_maxima();
  ivec3 hill_climb(ivec3 in);
  void highlight(std::vector<glm::ivec3> points);
  std::vector<ScaleBlob*> find_blobs();
  void draw_blobs(std::vector<ScaleBlob*>, const char *mode="gm");
  void add_blobs(std::vector<ScaleBlob*>);
  void color_blobs(std::vector<ScaleBlob*>, float color);
  ScaleBlob* compute_blob_tree();
  BSPTree<ScaleBlob> get_bsp(int depth);

  void capture(Nrrd *nin);
  void init(Nrrd *nin);
  Nrrd *commit(Nrrd *nout = 0);
  void destroy();

  static void print_kernel(DiscreteKernel k);
};

#endif