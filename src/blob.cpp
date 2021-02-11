#include "blob.h"
#include <cstdio>
#include <cmath>
#include <limits>
#include <LBFGS.h>
#include <iostream>
extern "C"{
  #include <lbfgs.h>
}
#define pi (3.14159265358979323846264338)

ScaleBlob::ScaleBlob(){

  position = vec3(0);
  
  shape = mat3x3(0);
  // eigs  = mat3x3(0);
  float inf = std::numeric_limits<float>::infinity();
  min   = vec3(inf, inf, inf);
  max   = vec3(0);

  parent = 0;
  scale  = -1;
  // volume = 0;

  npixels = 0;
  n      = 0;
  npass  = 0;

  peakvalue = 1;

  pred = std::vector<ScaleBlob*>();
  succ = std::vector<ScaleBlob*>();

  children = std::vector<ScaleBlob*>();

  model.type = '_';



  pass_sum_imagev = 0;
  pass_sum_expect = 0;
  // initialized = false;
}

static void printm(glm::mat3 m){
  printf("matrix"
    "\n   %.2f %.2f %.2f;"
    "\n   %.3f %.3f %.3f;"
    "\n   %.3f %.3f %.3f;"
    "\n",
    m[0][0],m[0][1],m[0][2],
    m[1][0],m[1][1],m[1][2],
    m[2][0],m[2][1],m[2][2]);
}

void ScaleBlob::pass(vec3 point, float value){
  if(npass == 0){  // calculate mean, min, max.
    npixels += 1;
    position += dvec3(point*value);
    n += value;
    if(point.x<min.x)min.x = point.x;
    if(point.y<min.y)min.y = point.y;
    if(point.z<min.z)min.z = point.z;
    if(point.x>max.x)max.x = point.x;
    if(point.y>max.y)max.y = point.y;
    if(point.z>max.z)max.z = point.z; 
  }
  if(npass == 1){  // calculate covariance.
    // printf("yo!");
    vec3 v = point - vec3(position);
    shape[0][0] += v.x*v.x*value/(n);
    shape[0][1] += v.x*v.y*value/(n);
    shape[0][2] += v.x*v.z*value/(n);
    shape[1][1] += v.y*v.y*value/(n);
    shape[1][2] += v.y*v.z*value/(n);
    shape[2][2] += v.z*v.z*value/(n);

  }
  if(npass == 2){
    // vec3 p = point - vec3(position);
    // pass_sum_expect += exp(-0.5 * glm::dot(p, invCov*p));
    // pass_sum_imagev += value;
    // // pass_distances.push_back();
    // pass_distances_weights.push_back(std::make_pair(glm::dot(p, invCov*p), value/n));
    // printf("%.2f; ", glm::dot(p, invCov*p));
  }
}
void ScaleBlob::commit(){
  // printf("commit %p\n", this);
  if(npass == 0){  // compute mean.
    position /= double(n);
    npass = 1;
  }
  else if(npass == 1){  // compute covariance matrix.
    // shape /= double(n);
    shape[1][0] = shape[0][1];
    shape[2][0] = shape[0][2];
    shape[2][1] = shape[1][2];
    // npass = 2;

    invCov    = glm::inverse(mat3(shape));
    detCov    = fabs(glm::determinant(shape));
    // // pdfCoef   = pow(glm::determinant(invCov*pi*2.0),-0.5);
    pdfCoef   = pow(2.0*pi, -1.5)*pow(glm::determinant(invCov), -0.5);

    // fshape = shape;

    covariance << shape[0][0], shape[0][1], shape[0][2],
                  shape[1][0], shape[1][1], shape[1][2],
                  shape[2][0], shape[2][1], shape[2][2];

    npass = 2;
    ll_is_blob = 0;
    GMM.gaus_scale = 1.0f/(pow(2*pi,3.0/2) * sqrt(detCov));
    printf("commit %f, %f\n", n, GMM.gaus_scale);
    printf("mean = %.4f %.4f %.4f\n", position.x, position.y, position.z);
    printm(shape);
  }
  else if(npass == 2){
    // float distsq95 = 1;
    // if(pass_distances_weights.size() > 0){

    //   std::sort(pass_distances_weights.begin(), pass_distances_weights.end());
    //   float sum = 0;
    //   int i=0;
    //   while(sum < 0.5){
    //     sum += pass_distances_weights[i].second;
    //     ++i;
    //   }
    //   int ind95 = i-1;
    //   // std::nth_element(pass_distances.begin(), pass_distances.begin() + ind95, pass_distances.end());
    //   float avg = 0;
    //   float qsum = 0;
    //   for(int i=0;i<pass_distances_weights.size(); i++){
    //     avg += (pass_distances_weights[i].first) * pass_distances_weights[i].second;
    //     if(pass_distances_weights[i].first < 1){
    //       qsum += pass_distances_weights[i].second;
    //     }
    //   }
    //   // avg /=   pass_distances.size();
    //   printf("sum = %.2f\n", qsum);
    //   // printf("n = %.3f\n", n);
    //   printf("ind95 %d/%d\n", ind95, pass_distances_weights.size());
    //   printf("avg = %.2f\n", avg);
    //   // printf("median = %.4f", 3.f*pow(1.-(2./(9.*3.  )), 3.f));
    //   distsq95 = pass_distances_weights[ind95].first;
    // }

    // printf("d = %.4f\n", sqrt(distsq95));

    // shape *= sqrt(distsq95);
    // min -= vec3(10,10,10);
    // max += vec3(10,10,10);

    // invCov    = glm::inverse(mat3(shape));
    // detCov    = fabs(glm::determinant(shape));
    // pdfCoef   = pow(2.0*pi, -1.5)*pow(glm::determinant(invCov), -0.5);

    // // // fshape = shape;

    // covariance << invCov[0][0], invCov[0][1], invCov[0][2],
    //               invCov[1][0], invCov[1][1], invCov[1][2],
    //               invCov[2][0], invCov[2][1], invCov[2][2];

    peakvalue = pass_sum_imagev / pass_sum_expect;
  // initialized = true;
  }
}
float ScaleBlob::pdf(vec3 p){
  // float vcenter = exp(-0.5 * 1);
  // float mul     = 1.f/exp(-0.5 * 1);
  // static float pdfCoef = n*pow(2.0*pi, -1.5)*pow(glm::determinant(invCov), -0.5);
  p = p - vec3(position);
  // float v = pdfCoef * exp(-0.5 * glm::dot(p,(invCov*p)));

  // the center should be 1.
  float v = peakvalue * exp(-0.5 * glm::dot(p, invCov*p));
  // printf("%.2f ", v );
  // if(v>1)v=1;
  return v;
}
float ScaleBlob::gauspdf(vec3 p){
  p = p - vec3(position);
  float pMp = glm::dot(p, invCov*p);
  if(pMp > 10)return 0;
  // printf("pmp = %f %f\n", pMp, GMM.gaus_scale);
      // printm(shape);
  float v = GMM.gaus_scale * exp(-0.5 * pMp);
  return v;
}
float ScaleBlob::gmmpdf(vec3 p){
  p = p - vec3(position);
  float pMp = glm::dot(p, invCov*p);
  if(pMp > 10)return 0;
  // printf("pmp = %f %f\n", pMp, GMM.gaus_scale);
      // printm(shape);
  float v = GMM.alpha * GMM.gaus_scale * exp(-0.5 * pMp);
  return v;
}
float ScaleBlob::cellpdf(vec3 p){
  p = p - vec3(position);
  float mag = glm::dot(p,(invCov*p));
  if(mag<0.1f)  return 1.f;
  if(mag<1.f)   return 0.8f;
  else          return 0.8f*(1 - 0.4 * (-1 + mag));
  // else         return float(erf(2-mag)*0.5+0.5);
  // float mag = glm::dot(p,(invCov*p));
  // return 1.f/(0.1f + 0.05f*mag*mag);
}
float ScaleBlob::ellipsepdf(vec3 p){
  p = p - vec3(position);
  float mag = glm::dot(p,(invCov*p));
  if(mag<0.1f)  return 1.f;
  if(mag<2.f)   return 0.8f;
  else          return 0.f;
  // else         return float(erf(2-mag)*0.5+0.5);
  // float mag = glm::dot(p,(invCov*p));
  // return 1.f/(0.1f + 0.05f*mag*mag);
}
float ScaleBlob::cellerf(vec3 p){
  p = p - vec3(position);
  float mag = glm::dot(p,(invCov*p));
  if(mag<1.f)  return 1.f;
  if(mag>3.5f) return 0.f;
  // else         return 1 - 0.4 * (-1 + mag);
  else         return float(erf(2-mag)*0.5+0.5);
  // float mag = glm::dot(p,(invCov*p));
  // return 1.f/(0.1f + 0.05f*mag*mag);
}
float ScaleBlob::celldot(vec3 p){
  p = p - vec3(position);
  float mag = sqrt(glm::dot(p,p));
  // if(mag < 1.f)printf("mag=%.2f\n", mag);
  if(mag<=1.42f)  return 1.f;
  if(mag>3.5f) return 0.f;
  else         return 0.f;
  // else         return 1 - 0.4 * (-1 + mag);
  // else         return float(erf(2-mag)*0.5+0.5);
  // float mag = glm::dot(p,(invCov*p));
  // return 1.f/(0.1f + 0.05f*mag*mag);
}
float ScaleBlob::outlinepdf(vec3 p){
  p = p - vec3(position);
  float mag = glm::dot(p,(invCov*p));
  if(mag<2.f)  return 0.f;
  if(mag<3.5f) return 0.5f;
  else         return 0;
  // if(mag>3.5f) return 0.f;
  // else         return 1 - 0.4 * (-1 + mag);
  // else         return float(erf(2-mag)*0.5+0.5);
  // float mag = glm::dot(p,(invCov*p));
  // return 1.f/(0.1f + 0.05f*mag*mag);
}
float ScaleBlob::erf_pdf(vec3 p){
  // model.alpha = 1;
  // model.beta = 4;
  // model.kappa = [0,1];

  // p = p-vec3(position);
  // float rr = glm::dot(p, (invCov*p));
  // vec3  gv = (invCov + glm::transpose(invCov)) * p;
  // float ss = (rr-1)/(0.000001 + glm::length(gv));
  // return (model.kappa * erf(-ss * model.alpha));
  // return 0;

  p = p - vec3(position);
  float mag = sqrt(glm::dot(p,(invCov*p)));
  if(mag<model.alpha){
    return model.kappa;
    // return model.kappa * 1.f;
  }
  else{
    mag -= model.alpha;
    mag *= model.beta;
    if(mag > 4){
      return 0;
    }
    return model.kappa * float(erf(2-mag)*0.5+0.5);
  }
  return 0;
}
float ScaleBlob::generalized_multivariate_gaussian_pdf(vec3 p){
  // if(model.kappa<0.1f)return celldot(p);
  // return 1.f;
  // printf("%.2f %.2f %.2f\n", alpha, beta, kappa);
  p = p-vec3(position);
  float mag = glm::dot(p,(invCov*p));
  float v = model.kappa*exp(-pow(model.alpha*mag, model.beta));
  // float v = exp(-0.5 * glm::dot(p, invCov*p));
  return v;
}
float ScaleBlob::modelpdf(vec3 p){
  if(model.type == 'e'){
    return erf_pdf(p);
  }else{
    return generalized_multivariate_gaussian_pdf(p);
  }
}
void ScaleBlob::print(){
  printf("blob %.2f at %.2f %.2f %.2f; xyz %.3f %.3f %.3f; xy/xz/yz %.3f %.3f %.3f\n", n,
    position[0],position[1],position[2],
    shape[0][0],shape[1][1],shape[2][2],
    shape[0][1],shape[0][2],shape[1][2]);
}
void ScaleBlob::printtree(int depth){
  // printf("\n");
  // for(int i=0;i<depth;++i)printf(" ");
  printf("(%.f",scale);
  for(auto c : children){
    c->printtree(depth+1);
  }
  printf(")");
} 
// float ScaleBlob::distance(vec3 p){
//   // = distance between (p, shape) and position
//   // = distance between (c, shape) and origin.
//   // vec3 c  = vec3(position) - p;
//   // vec3 x  = c + shape*vec3(0,0,1); // inital point x on ellipsoid.
//   // vec3 a(1,0,0);
//   return 0
// }

using namespace Eigen;

/** 
 * return xTAx
 */
static inline float ell(const Matrix3f &A, const Vector3f &x){
  float c00 = x[0]*x[0]*A(0,0);
  float c01 = x[0]*x[1]*A(0,1);
  float c02 = x[0]*x[2]*A(0,2);
  float c11 = x[1]*x[1]*A(1,1);
  float c12 = x[1]*x[2]*A(1,2);
  float c22 = x[2]*x[2]*A(2,2);

  return c00 + c11 + c22 + 2*c01 + 2*c02 + 2*c12;
}
float ScaleBlob::distance(ScaleBlob *blob){
  // printf("d");
  // printf("d...");
  /*** intersection distance between ellipsoids ***/
  using namespace Eigen;

  // define objective function.
  class EllipsoidDistanceObjective{
  public:
    EllipsoidDistanceObjective(const Matrix3f &A, const Vector3f p): A(A), p(p) {}
    Matrix3f A;
    Vector3f p;
    float operator()(const VectorXf &x, VectorXf &grad){
      float xTAx = ell(A, x);
      float xAx = ell(A, x);
      float xTx  = x.dot(x);
      float xp  = x.dot(p);
      float VxAx = sqrt(xAx);
      Vector3f Ax  = A*x;
      Vector3f xTAAT = (x.transpose() * (A + A.transpose())).transpose();

      grad = (2.f*x/(xAx)) + (-xTx/(xAx*xAx) * xTAAT) - (((1.f/VxAx)*2*p) + (2*xp)*(-xTAAT/(2.f*(VxAx*VxAx*VxAx))));
      // printf("] VxAx   = %.30f\n", VxAx);
      // printf("] p      = %.30f %.30f %.30f\n", p[0], p[1], p[2]);
      // printf("] 2/VxAx = %.30f\n", (1.f/VxAx)*2.f);
      // std::cout << "compute:\n"
      //   << (2.f*x/(xAx)) << "\n"
      //   << (-xTx/(xAx*xAx) * xTAAT) << "\n"
      //   << ((1.f/VxAx)*2.f*p) << "\n"
      //   << ((2.f*xp)*(-xTAAT/(2.f*(VxAx*VxAx*VxAx)))) << "\n";
      return (1.f/xAx * xTx) - (1.f/sqrt(xAx))*2.f*xp;
    }
    // map R^3 -> surface of ellipsoid.
    static inline Vector3f g(Matrix3f &A, const VectorXf &x){
      return x*(1.f/sqrt(ell(A,x)));
    }
  };

  // perform precomputation.
  typedef SelfAdjointEigenSolver<Eigen::Matrix3f> Solver3f;
  Vector3f p0(0, 0, 0);


  Solver3f s_A0(covariance);
  Solver3f s_A1(blob->covariance);

  Matrix3f VA0 = s_A0.operatorSqrt();
  Matrix3f VA1 = s_A1.operatorSqrt();

  Matrix3f VA0i = VA0.inverse();
  Matrix3f VA1i = VA1.inverse();

  Matrix3f VA0A1i  = VA0 * VA1i;
  Matrix3f VA1A0iv = VA1 * VA0i;
  
  // std::cout << "VA0 = \n" << VA0 << "\n\n";
  // std::cout << "VA1 = \n" << VA1 << "\n\n";

  // std::cout << "VA0i = \n" << VA0i << "\n\n";
  // std::cout << "VA1i = \n" << VA1i << "\n\n";

  // std::cout << "VA0A1i = \n" << VA0A1i  << "\n\n";
  // std::cout << "VA1A0i = \n" << VA1A0iv << "\n\n";

  Vector3f p1 = VA0 * Vector3f( position[0]-blob->position[0],
                                position[1]-blob->position[1],
                                position[2]-blob->position[2]) ;

  VectorXf x  = EllipsoidDistanceObjective::g(VA0A1i, p1);
  
  Matrix3f Aq = VA1A0iv * VA1A0iv;
  Matrix3f As;
  As <<   Aq(0,0),               (Aq(1,0)+Aq(0,1))/2.f, (Aq(2,0)+Aq(0,2))/2.f,
          (Aq(1,0)+Aq(0,1))/2.f, Aq(1,1),               (Aq(2,1)+Aq(1,2))/2.f,
          (Aq(2,0)+Aq(0,2))/2.f, (Aq(2,1)+Aq(1,2))/2.f, Aq(2,2);

  // printf("\n");
  // printf("p0 = %.3f, %.3f, %.3f\n", position.x, position.y, position.z);
  // printf("p1 = %.3f, %.3f, %.3f\n", blob->position.x, blob->position.y, blob->position.z);
  // std::cout << "A0   = \n" << covariance       << "\n\n";
  // std::cout << "A1   = \n" << blob->covariance << "\n\n";


  // std::cout << "p1   = \n" << p1 << "\n\n";
  // std::cout << "x    = \n" << x  << "\n\n";
  // std::cout << "Aq   = \n" << Aq  << "\n\n";
  // std::cout << "As   = \n" << As  << "\n\n";

  LBFGSpp::LBFGSParam<float> param;
  param.past  = 2;
  param.delta = 0.00001f;
  LBFGSpp::LBFGSSolver<float> solver(param);

  float fx;

  // printf("minimizing %f %f.\n", n, blob->n);
  EllipsoidDistanceObjective obj(As, p1);
  // printf("v");
  int nitr = solver.minimize(obj, x, fx);

  x = EllipsoidDistanceObjective::g(As, x);
  float distance = (x - p1).norm();
  // if(distance < 1.f){
  //   // printf("(%.3f %.3f %.3f) - (%.3f %.3f %.3f)\nd=%.2f\n", position[0], position[1], position[2], blob->position[0], blob->position[1], blob->position[2], distance);    
  // }
  // printf(".../d\n");
  return distance;  


  // compute distance between ellipsoids using newton's method.

  // vec3 p;


  /*** My distance: ***/
  // return (glm::length(blob->position - position)) + fabs(cbrt(detCov) - cbrt(blob->detCov));

  /*** Wasserstein metric:  ***/

  // using namespace Eigen;
  // SelfAdjointEigenSolver<Eigen::Matrix3f> solver(covariance);
  
  // Matrix3f sqc2 = solver.operatorSqrt();
  // Matrix3f c2c1c2 = sqc2 * blob->covariance * sqc2;
  
  // solver = SelfAdjointEigenSolver<Matrix3f>(c2c1c2);
  // Matrix3f sqrtc2c1c2 = solver.operatorSqrt();
  // Matrix3f whole = blob->covariance + covariance - (2.f * sqrtc2c1c2);

  // float trace = whole.trace();
  // vec3  delta = blob->position - position;
  // return (dot(delta, delta)) + trace;


  /*** simple distance ***/
  // return length(blob->position - position);
}

float ScaleBlob::wasserstein_distance(ScaleBlob* blob){
    /*** Wasserstein metric:  ***/

  using namespace Eigen;
  SelfAdjointEigenSolver<Eigen::Matrix3f> solver(covariance);

  // std::cout << covariance << std::endl << std::endl;
  
  Matrix3f sqc2 = solver.operatorSqrt();
  Matrix3f c2c1c2 = sqc2 * blob->covariance * sqc2;
  
  solver = SelfAdjointEigenSolver<Matrix3f>(c2c1c2);
  Matrix3f sqrtc2c1c2 = solver.operatorSqrt();
  Matrix3f whole = blob->covariance + covariance - (2.f * sqrtc2c1c2);

  float trace = whole.trace();
  vec3  delta = blob->position - position;
  return (dot(delta, delta)) + trace;
}

float ScaleBlob::covmaxev(){
  using namespace Eigen;
  SelfAdjointEigenSolver<Eigen::Matrix3f> solver(covariance);
  auto evs = solver.eigenvalues();
  // printf("evs %.2f %.2f %.2f %.2f\n", 1/evs[0], 1/evs[1], 1/evs[2], cbrt(detCov));
  // if(evs[0] != evs[0])return 0;
  float ev = 1.f/fabs(evs[0]);
  if(isnan(ev) || isinf(ev))return 0;
  return ev;
  // return 3.f*cbrt(detCov);
}

float ScaleBlob::covconditionnumber(){
  using namespace Eigen;
  SelfAdjointEigenSolver<Eigen::Matrix3f> solver(covariance);
  auto evs = solver.eigenvalues();
  // printf("evs %.2f %.2f %.2f %.2f\n", 1/evs[0], 1/evs[1], 1/evs[2], cbrt(detCov));
  // if(evs[0] != evs[0])return 0;
  float ev_max = 1.f/fabs(evs[0]);
  float ev_min = 1.f/fabs(evs[2]);
  if(isnan(ev_max) || isinf(ev_max))return 0;
  if(isnan(ev_min) || isinf(ev_min))return 0;
  return ev_max/ev_min;
  // return 3.f*cbrt(detCov);
}