#include <teem/meet.h>
#include <teem/nrrd.h>
#include <glm/glm.hpp>
#include "synth.h"
#include <vector>

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "filter.h"
#include "blob.h"

#define pi (3.14159265358979)
#define gaus(x,stdv) (exp(-((x)*(x))/(2*(stdv)*(stdv)))/((stdv)*sqrt(2*pi*pi*pi*2*2)))
#define ngaus(x,stdv) (exp(-((x)*(x))/(2*(stdv)*(stdv))))
// #define mgaus(x,stdv) (exp(-((x)*(x))))

#define DIM 100
#define SDIM 80
#define HDIM 200

// problems with >=85 and 10/20

void synth_ex(ArExperiment *e){
  using namespace glm;
  int a0=SDIM,a1=SDIM,a2=SDIM;
  short *data = new short[a0*a1*a2];

  ArFilter filter;
  float *fdata = new float[a0*a1*a2];
  float *fdata2 = new float[a0*a1*a2];

  Nrrd *nval = nrrdNew();
  Nrrd *ntemp = nrrdNew();

  srand ((unsigned long long)nval);
  // nrrdWrap_va(nval, fdata, nrrdTypeFloat, 3, a0,a1,a2);
  nrrdWrap_va(ntemp, fdata, nrrdTypeFloat, 3, a0,a1,a2);


  filter.init(ntemp);


  srand ((unsigned long long)nval);


  float sigma = 5.f;
  // DiscreteKernel blur = filter.gaussian(1.f, 5);
  for (int ti = 0; ti < (e->high - e->low) + 1; ti++){


    double x,y,z,v;

    // sphere.
    for(int i=0;i<a0*a1*a2;i++){
      x = (i%(SDIM));
      y = ((i/SDIM)%SDIM);
      z = (i/(SDIM*SDIM));

      vec3 p(x,y,z);
      vec3 q1 = p - vec3(SDIM/2,SDIM/2,SDIM/2) - vec3(0,0,0);
      // vec3 q2 = p - vec3(SDIM/2,SDIM/2,SDIM/2) + vec3(2,0,0);
      v = 0;

      v += ngaus((q1.x*q1.x + q1.y*q1.y + q1.z*q1.z), 14);
      // v += ngaus((q2.x*q2.x + q2.y*q2.y + q2.z*q2.z), 14);
      // double d = (q.x*q.x + q.y*q.y + q.z*q.z);
      // v = (d<100)?1:0;

      fdata[i] = v;
    }
    // filter.capture(ntemp);
  
    // DiscreteKernel kernel = filter.gaussian(sigma, int(sigma*3));
    // filter.set_kernel(kernel);
    // filter.filter();
    // filter.laplacian3d();
    // filter.normalize();
    // std::vector<glm::ivec3> maxima = filter.find_maxima();
    // printf("maxima = %d\n", maxima.size());
    // filter.highlight(maxima);

    // filter.commit(ntemp);

    // kernel.destroy();

    // convert float to short
    for(int i=0;i<a0*a1*a2;i++){
      data[i] = fdata[i] * 30000;
    }

    nrrdWrap_va(nval, data, nrrdTypeShort, 3, a0,a1,a2);

    printf("save to %s\n", e->getfilepath(e->low + ti).c_str());
    nrrdSave(e->getfilepath(e->low + ti).c_str(), nval, NULL);

    sigma += 0.5f;
  }


  filter.destroy();

  nrrdNix(nval);
  nrrdNix(ntemp);


  delete[] data;
  delete[] fdata;
}

void synth_highschool(ArExperiment *e){
  using namespace glm;
  int a0=HDIM,a1=HDIM,a2=HDIM;
  short *data = new short[a0*a1*a2];

  ArFilter filter;
  float *fdata = new float[a0*a1*a2];
  float *fdata2 = new float[a0*a1*a2];

  Nrrd *nval = nrrdNew();
  Nrrd *ntemp = nrrdNew();

  srand ((unsigned long long)nval);
  // nrrdWrap_va(nval, fdata, nrrdTypeFloat, 3, a0,a1,a2);
  nrrdWrap_va(ntemp, fdata, nrrdTypeFloat, 3, a0,a1,a2);
  DiscreteKernel kernel = filter.gaussian(1.f, 5);


  filter.init(ntemp);
  filter.set_kernel(kernel);


  srand ((unsigned long long)nval);



  for (int ti = 0; ti < 1; ti++){
    for(int i=0;i<a0*a1*a2; i++){
      fdata[i] = 0;
    }
    // filter.capture(ntemp);
    // filter.normalize();
    // filter.commit(ntemp);
    // printf("data: %p %p\n", fdata, ntemp->data);
    // memcpy(fdata, ntemp->data, sizeof(float) * a0*a1*a2);

    double x,y,z,v;
    float minv = 1000000;
    for(int i=0;i<a0*a1*a2;i++){
      x = (i%(HDIM)) - (HDIM/2);
      y = ((i/HDIM)%HDIM) - (HDIM/2);
      z = (i/(HDIM*HDIM)) - (HDIM/2);

      x *= 1.f;
      y *= 1.f;
      z *= 1.f;
      vec3 p(x,y,z);
      v = 0;

      vec3 v1(-2,  1, 2);
      vec3 v2(-2, -1, 2);
      vec3 v3( 1.5,  .8, 2);
      
      v1 = glm::normalize(v1);
      v2 = glm::normalize(v2);
      v3 = glm::normalize(v3);

      float t1=x, t2=y, t3=z;

      // v =  pow(glm::dot(t3-v3 - t1*v1, t2*v2 - t1*v1), 2);
      v += fabs(glm::dot((t3*v3) - (t1*v1), (t2*v2) - (t1*v1)));          // orthogonal
      v += fabs(sqrt(glm::dot(t1*v1-t2*v2, t1*v1-t2*v2))-10);     // v1, v2 = 10
      v += fabs(sqrt(glm::dot(t1*v1-t3*v3, t1*v1-t3*v3))-20);     // v1, v3 = 20
      v += fabs(sqrt(glm::dot(t1*v1-t3*v3, t1*v1-t3*v3))-22.36);     // v2, v3 = 22.36
      // v += glm::dot(t1*v1 -t2*v2, t1*v1 -t2*v2) - 10;
      // if(v<10)printf("%.2f; ", v);
      if(t1<0 || t2<0 || t3<0)v = 100;
      if(v < minv){
        printf("v= %.2f; %.2f %.2f %.2f\n", v, x, y, z);
        minv = v;
      }
      if(v<=20){
        v = 1.f/v;
      }else{
        v=0;
      }
      // if(v <= 1)v = 0;
      fdata[i] = v;
    }

    filter.capture(ntemp);
    filter.normalize();
    filter.commit(ntemp);

    // convert float to short
    for(int i=0;i<a0*a1*a2;i++){
      data[i] = fdata[i] * 30000;
    }

    printf("dims %d %d %d\n", a0, a1, a2);
    nrrdWrap_va(nval, data, nrrdTypeShort, 3, a0,a1,a2);
    // for(int i=e->low; i<=e->high; i++){
    printf("save to %s\n", e->getfilepath(e->low + ti).c_str());
    nrrdSave(e->getfilepath(e->low + ti).c_str(), nval, NULL);
    // }
  }


  kernel.destroy();
  filter.destroy();

  // nrrdSave("/home/ashwin/data2/synth/001.nrrd", nval, NULL);
  // nrrdSave("/home/ashwin/data2/synth/002.nrrd", nval, NULL);
  // nrrdSave("/home/ashwin/data2/synth/003.nrrd", nval, NULL);
  nrrdNix(nval);
  nrrdNix(ntemp);


  delete[] data;
  delete[] fdata;
}


void synth_single_varying(ArExperiment *e){
  using namespace glm;
  int a0=SDIM,a1=SDIM,a2=SDIM;
  short *data = new short[a0*a1*a2];

  ArFilter filter;
  float *fdata = new float[a0*a1*a2];
  float *fdata2 = new float[a0*a1*a2];

  Nrrd *nval = nrrdNew();
  Nrrd *ntemp = nrrdNew();

  srand ((unsigned long long)nval);
  // nrrdWrap_va(nval, fdata, nrrdTypeFloat, 3, a0,a1,a2);
  nrrdWrap_va(ntemp, fdata, nrrdTypeFloat, 3, a0,a1,a2);
  DiscreteKernel kernel = filter.gaussian(1.f, 5);


  filter.init(ntemp);
  filter.set_kernel(kernel);


  srand ((unsigned long long)nval);



  for (int ti = 0; ti < (e->high - e->low) + 1; ti++){



    // background noise
    for(int i=0;i<a0*a1*a2; i++){
      float v = float(rand() % 10000) / 10000.f;
      // if(v < 0.995) v=0;
      fdata[i] = v;
    }
    filter.capture(ntemp);
    // filter.filter();
    filter.normalize();
    filter.commit(ntemp);
    printf("data: %p %p\n", fdata, ntemp->data);
    // memcpy(fdata, ntemp->data, sizeof(float) * a0*a1*a2);


    std::vector<vec3> ps;
    // for(float xi=0.1f;xi<1;xi+=0.1f){
    {
      double x = SDIM * (0.3f + ti/15.f);
      double y = SDIM/2;
      double z = SDIM/2;
      printf("p %.3f %.3f %.3f\n",x,y,z);
      ps.push_back(vec3(x,y,z));
      // ps.push_back(vec3(x,y+15,z));
    }
    // }

    double stdv = 3;

    double x,y,z,v;
    for(int i=0;i<a0*a1*a2;i++){
      x = (i%(SDIM));
      y = ((i/SDIM)%SDIM);
      z = (i/(SDIM*SDIM));
      vec3 p(x,y,z);
      v = 0;

      for(int i=0;i<ps.size();i++){
        vec3 q = p-ps[i];
        // q.z *= 1.f;
        // double d = (q.x*q.x + q.y*q.y + q.z*q.z*0.25f);
        // v += exp(-0.5 * d);
        // v += gaus(d,stdv);
        // v += gaus(d,0.8);
        // v = y;

        double d = (q.x*q.x + q.y*q.y + 0.5*q.x*q.z + q.z*q.z*0.25f) * 0.02;
        // double d = (q.x*q.x + q.y*q.y + q.z*q.z) * 0.1;
        // v += exp(-0.5 * d);
        if(d < 4)v=1;
        // if(fabs(x-40) < 2 && fabs(y-40) < 2 && fabs(z-40) < 8)v=1;
        else v=0;
      }

      // v /= (ps.size()*gaus(0,stdv));

      // v = 0;

      // v = 0;

      // printf("%.2f %.2f %.2f\n", p.x, p.y, p.z);

      // if(v >= 1){
      //   printf("v %.3f\n", v);
      // }

      // if(v<0)v=0;
      // if(v>1)v=1;

      // v *= 0.4;

      // add ramp

      // v += 0.3f * (x/SDIM);

      // add low-frequency gaussian

      // v += 0.1 * gaus(sqrt(x*x + y*y + z*z), SDIM) / gaus(0.f, SDIM);
      float alpha = 0.0f;
      fdata[i] = fdata[i]*alpha + v*(1-alpha);

      // data[i] = 0;
    }
    // printf("~~~~~ HI : %.4f\n", gaus(0,3));

    // box blur radius increasing with depth

    // for(int x=0; x<SDIM; x++){
    //   for(int y=0; y<SDIM; y++){
    //     for(int z=0; z<SDIM; z++){
    //       int i = x + y*SDIM + z*SDIM*SDIM;
    //       float v = 0;
    //       float n = 0;

    //       int r = 1 + (float(x)/float(SDIM)) * 10;
    //       for(int xx = x-r; xx <= x+r; xx++){
    //         for(int yy = y-r; yy <= y+r; yy++){
    //           for(int zz = z-r; zz <= z+r; zz++){
    //             if(xx>=0 && yy>=0 && zz>=0 && xx<SDIM && yy<SDIM && zz<SDIM){
    //               v += fdata[xx + yy*SDIM + zz*SDIM*SDIM];
    //               n += 1;
    //             }
    //           }
    //         }
    //       }
    //       fdata2[i] = v/n; 
    //     }
    //   }
    // }


    // convert float to short
    for(int i=0;i<a0*a1*a2;i++){
      data[i] = fdata[i] * 30000;
    }

    // for (int i=0;i<SDIM * SDIM * SDIM; i++){
    //   data[i] = 0;
    // }
    // data[0] = 0;=
    // data[1] = 30000;
    printf("dims %d %d %d\n", a0, a1, a2);
    nrrdWrap_va(nval, data, nrrdTypeShort, 3, a0,a1,a2);
    // for(int i=e->low; i<=e->high; i++){
    printf("save to %s\n", e->getfilepath(e->low + ti).c_str());
    nrrdSave(e->getfilepath(e->low + ti).c_str(), nval, NULL);
    // }
  }


  kernel.destroy();
  filter.destroy();

  // nrrdSave("/home/ashwin/data2/synth/001.nrrd", nval, NULL);
  // nrrdSave("/home/ashwin/data2/synth/002.nrrd", nval, NULL);
  // nrrdSave("/home/ashwin/data2/synth/003.nrrd", nval, NULL);
  nrrdNix(nval);
  nrrdNix(ntemp);


  delete[] data;
  delete[] fdata;
}
void synthtracking(ArExperiment *e){
  using namespace glm;
  int a0=SDIM,a1=SDIM,a2=SDIM;
  short *data = new short[a0*a1*a2];
  Nrrd *nval = nrrdNew();

  srand ((unsigned long long)nval);

  for (int ti = 0; ti < (e->high - e->low) + 1; ti++){
    std::vector<vec3> ps;
    for(float xi=0.1f;xi<1;xi+=0.1f){
      double x = xi*a0;
      double y = 20 + ti * 3;
      double z = 20;
      printf("p %.3f %.3f %.3f\n",x*a0,y*a1,z*a2);
      ps.push_back(vec3(x,y,z));
      ps.push_back(vec3(y,x,z+20));
      ps.push_back(vec3(x,SDIM-y,z+40));
      ps.push_back(vec3(SDIM-y,x,z+60));
    }

    double stdv = 3;

    double x,y,z,v;
    for(int i=0;i<a0*a1*a2;i++){
      x = (i%(SDIM));
      y = ((i/SDIM)%SDIM);
      z = (i/(SDIM*SDIM));
      vec3 p(x,y,z);
      v = 0;

      for(int i=0;i<ps.size();i++){
        vec3 q = p-ps[i];
        // q.z *= 1.f;
        double d = sqrt(q.x*q.x + q.y*q.y + q.z*q.z);
        v += gaus(d,stdv);
        // v += gaus(d,0.8);
        // v = y;
      }

      v /= (ps.size()*gaus(0,stdv));

      // v = 0;

      // v = 0;

      // printf("%.2f %.2f %.2f\n", p.x, p.y, p.z);

      if(v >= 1){
        printf("v %.3f\n", v);
      }

      if(v<0)v=0;
      if(v>1)v=1;

      // v = blob.erf_pdf(p);

      data[i] = v*30000;

      // data[i] = 0;
    }
    // for (int i=0;i<SDIM * SDIM * SDIM; i++){
    //   data[i] = 0;
    // }
    // data[0] = 0;=
    // data[1] = 30000;
    printf("dims %d %d %d\n", a0, a1, a2);
    nrrdWrap_va(nval, data, nrrdTypeShort, 3, a0,a1,a2);
    // for(int i=e->low; i<=e->high; i++){
    printf("save to %s\n", e->getfilepath(e->low + ti).c_str());
    nrrdSave(e->getfilepath(e->low + ti).c_str(), nval, NULL);
    // }
  }

  // nrrdSave("/home/ashwin/data2/synth/001.nrrd", nval, NULL);
  // nrrdSave("/home/ashwin/data2/synth/002.nrrd", nval, NULL);
  // nrrdSave("/home/ashwin/data2/synth/003.nrrd", nval, NULL);
  nrrdNix(nval);

  delete[] data;
}

void synth(){
  // return;
  // 10 20 60
  using namespace glm;
  int a0=DIM,a1=DIM,a2=DIM;
  short *data = new short[a0*a1*a2];
  Nrrd *nval = nrrdNew();

  srand ((unsigned long long)nval);


  std::vector<vec3> ps;
  for(int i=0;i<10;i++){
    double x = (rand() % a0/10)/(a0/10.0);
    double y = (rand() % a1/10)/(a1/10.0);
    double z = (rand() % a2/10)/(a2/10.0);
    printf("p %.3f %.3f %.3f\n",x*a0,y*a1,z*a2);
    ps.push_back(vec3(x,y,z));
  }

  ScaleBlob blob;
  blob.position = vec3(50,50,50);
  blob.model.alpha = 0.1f;
  blob.model.beta = 0.1;
  blob.model.kappa = 0.8;
  blob.invCov = glm::mat3(0.01f,0,0, 0,0.02f,0, 0,0,0.04f);

  double stdv = 0.125;

  double x,y,z,v;
  for(int i=0;i<a0*a1*a2;i++){
    x = (i%(DIM));
    y = ((i/DIM)%DIM);
    z = (i/(DIM*DIM));
    vec3 p(x,y,z);
    v = 0;

    for(int i=0;i<ps.size();i++){
      vec3 q = p-ps[i];
      // q.z *= 1.f;
      double d = sqrt(q.x*q.x + q.y*q.y + q.z*q.z);
      v += gaus(d,stdv);
      // v += gaus(d,0.8);
      // v = y;
    }

    v /= (ps.size()*gaus(0,stdv));

    // v = 0;

    // printf("%.2f %.2f %.2f\n", p.x, p.y, p.z);

    if(v >= 1){
      printf("v %.3f\n", v);
    }

    if(v<0)v=0;
    if(v>1)v=1;

    v = blob.erf_pdf(p);

    data[i] = v*30000;
  }
  data[0] = 0;
  data[1] = 30000;
  printf("dims %d %d %d\n", a0, a1, a2);
  nrrdWrap_va(nval, data, nrrdTypeShort, 3, a0,a1,a2);
  nrrdSave("/home/ashwin/data/synth/000.nrrd", nval, NULL);
  nrrdSave("/home/ashwin/data/synth/001.nrrd", nval, NULL);
  nrrdSave("/home/ashwin/data/synth/002.nrrd", nval, NULL);
  nrrdSave("/home/ashwin/data/synth/003.nrrd", nval, NULL);
  nrrdNix(nval);
  delete[] data;
}