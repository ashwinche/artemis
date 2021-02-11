#ifndef BLOB_H
#define BLOB_H

#include <glm/glm.hpp>
#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <algorithm> 
#include <vector>
using namespace glm;


class ScaleBlob{
public:
  ScaleBlob *parent;                // scale-tree parent
  std::vector<ScaleBlob*> children; // scale-tree children

  std::vector<ScaleBlob*> pred;     // temporal predecessors
  std::vector<ScaleBlob*> succ;     // temporal successors

  int    imode;      // the local maximum which seeded this blob.
  vec3    mode;       // the local maximum which seeded this blob.
  float peakvalue;  // value of the image at the peak.
  int npixels;      // number of pixels in blob.
  dvec3  position;  // mean of this blob in space.
  mat3x3 shape;    // covariance matrix of blob.
  int timestep;     // the timestep in which this blob exists.

  // temporary variables when constructing blob
  float pass_sum_imagev;
  float pass_sum_expect;
  float pass_sumsquaredweights;
  std::vector<std::pair<float,float>> pass_distances_weights;
  // mat3x3 fshape;    // covariance matrix of blob.

  struct{
    float nk;
    float alpha;
    float gaus_scale;
    std::vector<int>   pixels;
    std::vector<float> weights;
  }GMM;

  struct{
    float alpha;      // scaling parameter. 
    float beta;       // kurtosis parameter. 
    float kappa;      // magnitude paramater.
    char  type;       // either 'e'rf or 'g'aussian.
    float n;          // number of pixels explained by model.
    float error;      // total amount of error over these pixels.
                      // average error is error/n.
    vec3 min;         // min/max of portion of image explained 
    vec3 max;         // by model.
  }model;
  Eigen::Matrix3f covariance;


  // dmat3x3  eigs; // eigenvectors of covariance matrix.
  mat3x3 invCov;    // inverse of covariance matrix.
  float  detCov;    // determinant of covariance matrix.
  float  pdfCoef;   // |2*pi*covariance_matrix|^(-0.5)
  vec3   min;       // min bounding box of blob.
  vec3   max;       // max bounding box of blob.
  float scale;      // scale at which is blob was found.

  // double volume;
  float n;
  int npass;

  float ll_is_blob;  // log-likelihood that this blob is a gaussian.


  // bool initialized;

  ScaleBlob();
  float pdf(vec3 p);
  float gauspdf(vec3 p);
  float gmmpdf(vec3 p);
  float cellpdf(vec3 p);
  float cellerf(vec3 p);
  float celldot(vec3 p);
  float outlinepdf(vec3 p);
  float ellipsepdf(vec3 p);
  float generalized_multivariate_gaussian_pdf(vec3 p);
  float erf_pdf(vec3 p);

  float modelpdf(vec3 p);

  void pass(vec3 point, float value);
  void commit();
  void print();
  void printtree(int depth=0);
  float distance(ScaleBlob*);
  float wasserstein_distance(ScaleBlob*);

  float covmaxev();
  float covconditionnumber();
};

typedef ScaleBlob Blob;
#endif