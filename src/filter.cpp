#include "filter.h"
#include "HuangQS.h"
#include <vector>
#include <queue> 
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <unistd.h>
#include <climits>
#include <cmath>
#include <limits>
#include "util.h"
#include <random>
#include <time.h>
#define pi (3.14159265358979323846264)
#define gaus(x,stdv) (exp(-((x)*(x))/(2*(stdv)*(stdv)))/((stdv)*sqrt(2*pi)))
#define gausd2(x,sig)(exp(-x*x/(2.0*sig*sig))*(x*x-sig*sig) /       \
     (sig*sig*sig*sig*sig*2.50662827463100050241))


DiscreteKernel::DiscreteKernel(){
  // printf("data=0\n");
  data = 0;
  // printf("segfault\n");
  temp = 0;
}
DiscreteKernel::~DiscreteKernel(){
  // if(data) delete[] data;
  // if(temp) delete[] temp;
}
void DiscreteKernel::destroy(){
  if(data) delete[] data;
  if(temp) delete[] temp;
  data = 0;
  temp = 0;
}
struct thread_conv2d_info_shared{
  float *in;
  float *out;
  int xbegin;
  int xlen;
  int ylen;
  int zlen;
  int xstep;
  int ystep;
  int zstep;
  DiscreteKernel kernel;
};
struct tconv2d_data{
  float *in;
  float *out;

  int xstart;

  int xend;
  int ylen;
  int zlen;
  int xstep;
  int ystep;
  int zstep;
  DiscreteKernel kernel;
};
static void* tconv2d(void *tinfo){
  tconv2d_data *data = (tconv2d_data*) tinfo;
  
  float *in = data->in;
  float *out = data->out;
  
  int xstart = data->xstart;

  int xend = data->xend;
  int ylen = data->ylen;
  int zlen = data->zlen;
  int xstep = data->xstep;
  int ystep = data->ystep;
  int zstep = data->zstep;
  // printf("kernel = data->kernel.\n");
  DiscreteKernel kernel = data->kernel;
  // printf("success.\n");
  int skip  = zstep;
  int n = zlen;

  for(int x=xstart; x<xend; ++x){
    // printf("x=%d\n", x);
    for (int y=0; y<ylen; ++y){
      int start = x*xstep + y*ystep;
      double v = 0;
      
      // handle leading edge of input (within radius of edge).
      int i = 0;
      for(i=0;i<kernel.radius;++i){
        v=0;
        for(int j=0,ji=start+skip*(i-kernel.radius);j<kernel.support;++j,ji+=skip){
          if(ji<start){   // out of bounds. bleed.
            v += kernel.data[j] * in[start];
          }
          else if(ji<=start+skip*(n-1)){        // in bounds.
            v += kernel.data[j] * in[ji];
          }
        }
        out[start+skip*i] = float(v);
      }

      // handle trailing edge of input.
      for(i=n-kernel.radius;i<n;++i){
        v=0;
        for(int j=0,ji=start+skip*(i-kernel.radius);j<kernel.support;++j,ji+=skip){
          if(ji>start+skip*(n-1)){  // out of bounds. bleed.
            v += kernel.data[j] * in[start+skip*(n-1)];
          }
          else if(ji>=0){        // in bounds.
            v += kernel.data[j] * in[ji];
          }
        }
        out[start+skip*i] = float(v);
      }

      // convolve the center of the image.
      for(i=kernel.radius;i<n-kernel.radius;++i){
        v=0;
        for(int j=0,ji=start+skip*(i-kernel.radius);j<kernel.support;++j,ji+=skip){
          v += kernel.data[j] * in[ji];
        }
        out[start+skip*i] = float(v);
      }
    }
  }
  // printf("end thread.\n");
  return 0;
}
void ArFilter::conv2d(float *in, float *out, int xlen, int ylen, int zlen, int xstep, int ystep, int zstep, DiscreteKernel kernel){
  const int nthreads = 1;
  pthread_t threads[nthreads];
  tconv2d_data data[nthreads];

  for(int i=0;i<nthreads;++i){
    // printf("init thread %d\n", i);
    data[i].in = in;
    data[i].out = out;
    data[i].xstart = (xlen*i)/nthreads;
    data[i].xend = (xlen*(i+1))/nthreads;
    data[i].ylen = ylen;
    data[i].zlen = zlen;
    data[i].xstep = xstep;
    data[i].ystep = ystep;
    data[i].zstep = zstep;
    data[i].kernel = kernel;

    if(pthread_create(threads+i, NULL, tconv2d, data+i)){
      fprintf(stderr, "error creating conv2d thread.\n");
      exit(0);
    }
  }
  // printf("waiting...\n");
  for(int i=0;i<nthreads;++i){
    // printf("joining...\n");
    if(int err = pthread_join(threads[i], 0)){
      fprintf(stderr, "error joining conv2d thread. %d\n", err);
      exit(0);
    }
    // printf("joined.\n");
  }
  // printf("done.\n");
}
void ArFilter::filter(){
  // printf("filter.\n");
  // printf("curr=%d\n", self.curr);
  float *in  = self.buff[self.curr];
  float *out = self.buff[itempbuf()];

  // printf("conv: %p -> %p\n", in, out);
  conv2d(in, out, self.a0, self.a1, self.a2, self.w0, self.w1, self.w2, self.kernel);  // xy(z) axis.
  // printf("conv: %p -> %p\n", out, in);
  conv2d(out, in, self.a0, self.a2, self.a1, self.w0, self.w2, self.w1, self.kernel);  // x(y)z axis.
  // printf("conv: %p -> %p\n", in, out);
  conv2d(in, out, self.a1, self.a2, self.a0, self.w1,  self.w2, self.w0, self.kernel);  // (x)yz axis.

  int from,to;

  // for(int z=self.a2-1;z>=0;z--){
  //   for(int y=self.a1-1;y>=0;y--){
  //     for(int x=self.a0-1;x>=0;x--){
  //       to = x*self.w0 + y*self.w1 + z*self.w2;
        
  //       if(  x<self.kernel.radius || y<self.kernel.radius || z<self.kernel.radius
  //         || x>=self.a0-self.kernel.radius || y>=self.a1-self.kernel.radius || z>=self.a2-self.kernel.radius){
  //           out[to] = 0;
  //         continue;
  //       }
        
  //       from = (x-self.kernel.radius)*self.w0 + (y-self.kernel.radius)*self.w1 + (z-self.kernel.radius)*self.w2;
  //       out[to] = out[from];
  //     }
  //   }
  // }
  self.curr = itempbuf();
}
DiscreteKernel ArFilter::gaussian(double sigma, int radius, int d){
  // printf("create gaussian kernel: %.2f, %d\n",sigma, radius);
  if(sigma == 0){
    DiscreteKernel k;
    k.radius = 1;
    k.support = 3;
    k.data = new double[3];
    k.temp = new double[3];
    k.data[0] = 0;
    k.data[1] = 1;
    k.data[2] = 0;
    return k;
  }
  DiscreteKernel k;
  k.radius  = radius;
  k.support = radius*2 + 1;
  k.data = new double[k.support];
  k.temp = new double[k.support];
  int h = radius;
  if(d==2)  k.data[h] = gausd2(0,sigma);
  else      k.data[h] = gaus(0,sigma);
  for(int i=1;i<radius+1;++i){
    if(d==2)  k.data[h+i] = gausd2(i,sigma);
    else      k.data[h+i]   = gaus(i,sigma);
    k.data[h-i] = k.data[h+i];
  }
  // normalize so that the kernel sums to 1.0.
  double sum = 0;
  for(int i=0;i<k.support;++i)sum += k.data[i];
  for(int i=0;i<k.support;++i)k.data[i] /= sum;

  // print_kernel(k);
  return k;
}
DiscreteKernel ArFilter::interpolation(){
  DiscreteKernel k;
  k.radius = 1;
  k.support = 3;
  k.data = new double[3];
  k.temp = new double[3];
  k.data[0] = 0.5;
  k.data[1] = 0;
  k.data[2] = 0.5;
  return k;
}
void ArFilter::set_kernel(DiscreteKernel k){
  this->self.kernel = k;
}
void ArFilter::print_kernel(DiscreteKernel k){
  printf("kernel...\n  ");
  for(int i=0;i<k.support;++i){
    printf("%2.3f ",k.data[i]);
  }
  printf("\n/kernel...\n");
}
// struct filter_info{
//   int len;
//   float *in;
//   float *out;
// };
// static filter_info get_info(Nrrd *nin, Nrrd *nout, int channel){
//   filter_info fi;
//   NrrdAxisInfo *a = nin->axis;
//   fi.len = a[0].size * a[1].size * a[2].size *a[3].size;
//   fi.in = (float*)nin->data;
//   fi.out  = (float*)nout->data;
//   return fi;
// }
// void ArFilter::positive(Nrrd *nin, Nrrd *nout, int channel){
//   filter_info f = get_info(nin, nout, channel);
//   for(int i=channel;i<f.len;i+=2){
//     if(f.in[i]>0)f.out[i]=f.in[i];
//     else f.out[i] = 0;
//   }
// }
// void ArFilter::negative(Nrrd *nin, Nrrd *nout, int channel){
//   filter_info f = get_info(nin, nout, channel);
//   for(int i=channel;i<f.len;i+=2){
//     if(f.in[i]<0)f.out[i]= f.in[i];
//     else f.out[i] = 0;
//   }
// }
// void ArFilter::binary(Nrrd *nin, Nrrd *nout, int channel){
//   filter_info f = get_info(nin, nout, channel);
//   for(int i=channel;i<f.len;i+=2){
//     if(f.in[i]>0)f.out[i]=100;
//     else f.out[i] = 0;
//   }
// }

// double ArFilter::comp_max_laplacian(float *in){
//   double max = 0;
//   double lap;
//   int xi, xm;
//   for(int x=boundary;x<self.a0;++x){
//     for(int y=boundary;y<self.a1;++y){
//       xi = x*self.w0 + y*self.w1 + 1*self.w2;
//       for(; xi<self.w3-self.w2; xi+=self.w2){
//         lap = 2.0*in[xi] - in[xi-self.w2] - in[xi+self.w2];
//         if(lap>max)max=lap;
//       }
//     }
//   }
//   for(int x=boundary;x<self.a0;++x){
//     for(int z=boundary;z<self.a2;++z){
//       xi = x*self.w0 + 1*self.w1           + z*self.w2;
//       xm = x*self.w0 + (self.a1-1)*self.w1 + z*self.w2;
//       for(; xi<xm; xi+=self.w1){
//         lap = 2.0*in[xi] - in[xi-self.w1] - in[xi+self.w1];
//         if(lap>max)max=lap;
//       }
//     }
//   }

//   for(int y=boundary;y<self.a1;++y){
//     for(int z=boundary;z<self.a2;++z){
//       xi = 1*self.w0 + y*self.w1           + z*self.w2;
//       xm = (self.a0-1)*self.w0 + y*self.w1 + z*self.w2;
//       for(; xi<xm; xi+=self.w0){
//         lap = 2.0*in[xi] - in[xi-self.w0] - in[xi+self.w0];
//         if(lap>max)max=lap;
//       }
//     }
//   }
//   return 3.0*max/2.0;
// }


void ArFilter::deconvolve(){

  Deconvolve deconv;
  Imdim dim(self.a0, self.a1, self.a2);
  // PSF psf;
  float* I      = self.buff[self.curr];
  float* II     = self.buff[itempbuf()];
  float* orig   = new float[self.w3];
  for(int i=0;i<self.w3;i++)orig[i] = I[i];

  // for testing, blur I with the kernel. 
  // PSF psf(glm::mat3x3(0.4,0,0, 0,0.4,0, 0,0,0.4));
  // for(int x=0;x<self.a0;x++){
  //   for(int y=0;y<self.a1;y++){
  //     for(int z=0;z<self.a2;z++){
  //       int i = x + y*self.w1 + z*self.w2;
  //       II[i] = psf.convolve(I, x,y,z, dim);
  //     }
  //   }
  // }
  // self.curr = itempbuf();
  
  // I = II;


  float* Sn     = self.buff[itempbuf()];
  float* In     = new float[self.w3];
  float* C      = new float[self.w3];

  // for(int i=0;i<self.w3;i++)Sn[i] = 1.f;
  for(int i=0;i<self.w3;i++)Sn[i] = 1.f;

  printf("original image @ %p\n", I);
  printf("guess    image @ %p\n", Sn);
  printf("temp     image @ %p\n", C);


  for(int lritr = 0; lritr < 3; ++lritr){
    printf("itr %d\n", lritr);

    // In = Sn convolve P
    PSF psf(glm::mat3x3(0.2,0,0, 0,0.2,0, 0,0,0.2));
    for(int x=0;x<self.a0;x++){
      for(int y=0;y<self.a1;y++){
        for(int z=0;z<self.a2;z++){
          int i = x + y*self.w1 + z*self.w2;
          In[i] = psf.convolve(Sn, x,y,z, dim);
        }
      }
    }

    // In = I / In
    float error = 0;
    for(int i=0;i<self.w3;i++){
      In[i] = I[i] / In[i];
      error += In[i];
    }
    error /= self.w3;
    printf("error = %f\n", error);



    // C = In convolve P_reversed
    PSF psfT(-glm::mat3x3(0.2,0,0, 0,0.2,0, 0,0,0.2));
    for(int x=0;x<self.a0;x++){
      for(int y=0;y<self.a1;y++){
        for(int z=0;z<self.a2;z++){
          int i = x + y*self.w1 + z*self.w2;
          C[i] = psfT.convolve(In, x,y,z, dim);
          // C[i] = In[i];
          // printf("%.5f -> %.5f\n", C[i], In[i]);
        }
      }
    }

    // Sn = C * Sn
    for(int i=0;i<self.w3;i++)Sn[i] = C[i] * Sn[i];
  }

  // for(int i=0;i<self.w3/3;i++)Sn[i] *= 0.1;

  printf("sanity check.\n");

  // I = Sn concolve P
  PSF psf(glm::mat3x3(0.2,0,0, 0,0.2,0, 0,0,0.2));
  for(int x=0;x<self.a0;x++){
    for(int y=0;y<self.a1;y++){
      for(int z=0;z<self.a2;z++){
        int i = x + y*self.w1 + z*self.w2;
        I[i] = psf.convolve(Sn, x,y,z, dim);
      }
    }
  }


  delete[] In;
  delete[] C;
  delete[] orig;
  self.curr = itempbuf();
  normalize();
  printf("done.\n");
  // printf("Sn     @ %p\n", Sn);
  // printf("output @ %p\n", (self.buff[self.curr]));
}
void ArFilter::laplacianmasked(float scale){
  float *raw = self.buff[self.curr];
  int iraw = self.curr;
  laplacian3d();
  float *lap = self.buff[self.curr];
  for(int i=0;i<self.w3;i++){
    if(lap[i] == 0)raw[i] = 0;
  }
  self.curr = iraw;

}
struct tlaplacian3d_data{
  int zbegin;
  int zend;
  int boundary;
  ArFilter *filter;
  float *in;
  float *out;
};
static void* tlaplacian3d(void *vdata){
  tlaplacian3d_data *data = (tlaplacian3d_data*) vdata;
  
  ArFilter *filter = data->filter;
  float *in = data->in;
  float *out = data->out;
  int boundary = data->boundary;

  int xo;

  for(int z=data->zbegin;z<data->zend; z++){
    for(int y=0; y<filter->self.a1; y++){
      for(int x=0; x<filter->self.a0; x++){
        xo = x*filter->self.w0 + y*filter->self.w1 + z*filter->self.w2;

        if(x<1+boundary || y<1+boundary || z<1+boundary || x>=filter->self.a0-1-boundary || y>=filter->self.a1-1-boundary || z>=filter->self.a2-1-boundary){
          out[xo]=0;
          continue;
        }
        
        int p1=0,p2=0,p3=0,
            p4=0,p5=0,p6=0;
        double v;

        int xi = xo;
        p1 = xi-filter->self.w0;
        p2 = xi+filter->self.w0;
        p3 = xi-filter->self.w1;
        p4 = xi+filter->self.w1;
        p5 = xi-filter->self.w2;
        p6 = xi+filter->self.w2;

        // printf("%d %d -> %d %d\n",xi, in[xi], p1, in[p1]);

        double l1=0,l2=0,l3=0;

        l1 = 2*in[xi] - in[p1] - in[p2];
        l2 = 2*in[xi] - in[p3] - in[p4];
        l3 = 2*in[xi] - in[p5] - in[p6];
        // double mx = fmax(l1,fmax(l2,l3));
        // double mn = fmin(l1,fmin(l2,l3));
        // printf("%f %f %f -> %f %f\n",l1,l2,l3,mx,mn);
        v = (l1+l2+l3);
        // v *= 30000.0/max_laplacian;
        // v = fabs(in[p5]);
        if(v<0)v=0;
        out[xo] = float(v);
      }
    }
  }
  return 0;
}
void ArFilter::pass(){
  float *in  = self.buff[self.curr];
  float *out = self.buff[itempbuf()];
  memcpy(out, in, self.w3 * sizeof(float));
  self.curr = itempbuf();
}
void ArFilter::laplacian3d(int boundary){
  float *in  = self.buff[self.curr];
  float *out = self.buff[itempbuf()];

  

  const int nthreads = 24;
  pthread_t threads[nthreads];
  tlaplacian3d_data data[nthreads];

  for(int i=0;i<nthreads;++i){
    data[i].zbegin = (this->self.a2*i)/nthreads;
    data[i].zend = (this->self.a2*(i+1))/nthreads;
    data[i].boundary = boundary;
    data[i].filter = this;
    data[i].in = in;
    data[i].out = out;

    if(pthread_create(threads+i, NULL, tlaplacian3d, data+i)){
      fprintf(stderr, "error creating conv2d thread.\n");
      exit(0);
    }
  }
  for(int i=0;i<nthreads;++i){
    if(int err = pthread_join(threads[i], 0)){
      fprintf(stderr, "error joining conv2d thread. %d\n", err);
      exit(0);
    }
  }

  // printf("lap %p -> %p\n",in,out);
  // printf("nb %d %p %p\n",self.nbuf, self.buff[0], self.buff[1]);

  // printf("dims %d %d %d %d\n",self.a0,self.a0,self.a1,self.a2);
  // printf("max %d\n",self.w3);

  // double max_laplacian = comp_max_laplacian(in);

  self.curr = itempbuf();
}

void ArFilter::hessian3d(int boundary){
  float *in  = self.buff[self.curr];
  float *out = self.buff[itempbuf()];

  // printf("lap %p -> %p\n",in,out);
  // printf("nb %d %p %p\n",self.nbuf, self.buff[0], self.buff[1]);

  // printf("dims %d %d %d %d\n",self.a0,self.a0,self.a1,self.a2);
  // printf("max %d\n",self.w3);

  // double max_laplacian = comp_max_laplacian(in);
  int xo;
  for(int z=0;z<self.a2; z++){
    for(int y=0; y<self.a1; y++){
      for(int x=0; x<self.a0; x++){
        xo = x*self.w0 + y*self.w1 + z*self.w2;

        if(x<1+boundary || y<1+boundary || z<1+boundary || x>=self.a0-1-boundary || y>=self.a1-1-boundary || z>=self.a2-1-boundary){
          out[xo]=0;
          continue;
        }
        
        int p1=0,p2=0,p3=0,
            p4=0,p5=0,p6=0;
        double v;

        int xi = xo;

        int neighbor[27];
        int *n = neighbor;
        for(int i=0;i<27;i++){
          neighbor[i] = xi;
        }
        int nindex = 0;
        for(int zz=-1;zz<=1;++zz){
          for(int yy=-1;yy<=1;++yy){
            for(int xx=-1;xx<=1;++xx){
              int xii = xi + xx*self.w0 + yy*self.w1 + zz*self.w2;
              if(xii>=0 && xii < self.w3){
                neighbor[nindex] = xii;
              }
                ++nindex;
            }
          }
        }

        // printf("%d %d -> %d %d\n",xi, in[xi], p1, in[p1]);

        // double l1=0,l2=0,l3=0;

        float h11 = 2*in[xi] - in[neighbor[12]] - in[neighbor[14]];
        float h22 = 2*in[xi] - in[neighbor[10]] - in[neighbor[16]];
        float h33 = 2*in[xi] - in[neighbor[4]] - in[neighbor[22]];

        // 9, 11, 15, 17
        float h12 = in[n[15]] - in[n[9]] - in[n[17]] + in[n[11]];
        float h13 = in[n[3]] - in[n[5]] - in[n[21]] + in[n[23]];
        float h23 = in[n[1]] - in[n[7]] - in[n[19]] + in[n[25]];

        h11/=2;
        h22/=2;
        h33/=2;
        h12/=4;
        h13/=4;
        h23/=4;

        float a = h11;
        float b = h12;
        float c = -h13;
        float d = h22;
        float e = -h23;
        float f = h33;

        // a=0;
        // b=0;
        // c=0;
        // d=0;
        // e=0;
        // f=0;


        // (15 - 17) - (9 - 11) = 15 - 17 - 9 + 11
        // (5 - 3) - (23 - 21)
        // (7 - 1) - (25 - )

        // v = sqrt(a*a + b*b + c*c + d*d + e*e + f*f);

        v = a*(d*f - e*e) + b*(c*e - b*f) + c*(b*e - d*c);

        // a = a;
        // d = e;

        // h12 = h12/2;

        // glm::mat2x2 hessian(h11, h12, h12, h22);
        Eigen::Matrix<float, 3, 3> m;
        m << h11, h12, h13, 
             h12, h22, h23, 
             h13, h23, h33;
        
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float,3,3>> es(m);
        Eigen::Vector3f evs = es.eigenvalues();

        float ev1 = evs[0];
        float ev2 = evs[1];
        float ev3 = evs[2];
        // v = glm::determinant(hessian);

        a = h11;
        c = h12;
        b = h22;
        // float ev1 = 0.5*(a + b + sqrt((a-b)*(a-b) + 4*c*c));
        // float ev2 = 0.5*(a + b - sqrt((a-b)*(a-b) + 4*c*c));
        // v = h11*h22 - h12*h12;
        v = fmin(ev1, fmin(ev2, ev3));
        // v = fmax(ev1, fmax(ev2, ev3));
        // v = fmin(fabs(ev1), fmin(fabs(ev2), fabs(ev3)));
        // v = fmax(fabs(ev1), fmax(fabs(ev2), fabs(ev3)));

        // v = fabs(ev1*ev2*ev3);
        // if(ev1<-0.0001 || ev2<-0.0001 || ev3<-0.0001)v=0;
        // if(ev3<0)v=0;
        // v= fabs(ev1*ev2*ev3);
        v = fmax(v,0);
        // if(in[xi] > 0.8f && v > 0.01)printf("+ %.2f %.2f %.2f\n", ev1, ev2, ev3);
        // if(in[xi] < 0.2f && v > 0.01)printf("- %.2f %.2f %.2f\n", ev1, ev2, ev3);

        // v = ev1 * ev2 * ev3;
        // v = ev1 + ev2 + ev3;
        // v = h11 + h22;


        // v = -v;
        // v = max(h11, max(h22, h33));
        // double mx = fmax(l1,fmax(l2,l3));
        // double mn = fmin(l1,fmin(l2,l3));
        // printf("%f %f %f -> %f %f\n",l1,l2,l3,mx,mn);
        // v = (h12);
        // v *= 30000.0/max_laplacian;
        // v = fabs(in[p5]);
        // if(v<=0)v=0;
        // else v = in[xi];
        // else v = 1.f/v;
        // else v = ev1+ev2+ev3;
        out[xo] = float(v);
      }
    }
  }
  self.curr = itempbuf();
}

inline float fmedian(float a, float b, float c){
  if(a<b){            // a b
    if(b<c)return b;  // a b c
    if(c<a)return a;  // c a b
    else   return c;  // a c b
  }else{              // b a
    if(a<c)return a;  // b a c
    if(c<b)return b;  // c b a
    else   return c;  // b c a
  }
}
void ArFilter::max1(){
  float *in  = self.buff[self.curr];
  float *out = self.buff[itempbuf()];

  // printf("max1 %p -> %p\n", in, out);
  int xi,xm;

  // for(int x=0;x<self.a0;++x){
  //   for(int y=0;y<self.a1;++y){
  //     xi = x*self.w0 + y*self.w1 + 1*self.w2;
  //     for(; xi<self.w3-self.w2; xi+=self.w2){
  //       // out[xi] = 1.f;
  //       // out[xi+1] = in[xi];
  //       // out[xi] = in[xi];
  //       // out[xi] = in[xi+self.w2];
  //       out[xi] = fmedian(in[xi], in[xi-self.w2],in[xi+self.w2]);
  //       // out[xi] = max(in[xi],max(in[xi-self.w2],in[xi+self.w2]));
  //     }
  //   }
  // }
  // for(int x=0;x<self.a0;++x){
  //   for(int z=0;z<self.a2;++z){
  //     xi = x*self.w0 + 1*self.w1           + z*self.w2;
  //     xm = x*self.w0 + (self.a1-1)*self.w1 + z*self.w2;
  //     for(; xi<xm; xi+=self.w1){
  //       // in[xi] = out[xi];
  //       in[xi] = fmedian(out[xi],out[xi-self.w1],out[xi+self.w1]);
  //       // in[xi] = max(out[xi],max(out[xi-self.w1],out[xi+self.w1]));
  //     }
  //   }
  // }

  // for(int y=0;y<self.a1;++y){
  //   for(int z=0;z<self.a2;++z){
  //     xi = 1*self.w0 + y*self.w1           + z*self.w2;
  //     xm = (self.a0-1)*self.w0 + y*self.w1 + z*self.w2;
  //     for(; xi<xm; xi+=self.w0){
  //       // out[xi] = in[xi];
  //       out[xi] = fmedian(in[xi],in[xi-self.w0],in[xi+self.w0]);
  //       // out[xi] = max(in[xi],max(in[xi-self.w0],in[xi+self.w0]));
  //     }
  //   }
  // }
  for(int z=0;z<self.a2; z++){
    for(int y=0; y<self.a1; y++){
      for(int x=0; x<self.a0; x++){

        xi = x*self.w0 + y*self.w1 + z*self.w2;
        
        int neighbor[27];
        for(int i=0;i<27;i++){
          neighbor[i] = xi;
        }
        int nindex = 0;
        for(int xx=-1;xx<=1;++xx){
          for(int yy=-1;yy<=1;++yy){
            for(int zz=-1;zz<=1;++zz){
              int xii = xi + xx*self.w0 + yy*self.w1 + zz*self.w2;
              if(xii>=0 && xii < self.w3){
                neighbor[nindex] = xii;
                ++nindex;
              }
            }
          }
        }
        double v;
        v = fmedian(
          fmedian(
            fmedian(in[neighbor[0]], in[neighbor[1]], in[neighbor[2]]),
            fmedian(in[neighbor[3]], in[neighbor[4]], in[neighbor[5]]),
            fmedian(in[neighbor[6]], in[neighbor[7]], in[neighbor[8]])),
          fmedian(
            fmedian(in[neighbor[9]], in[neighbor[10]], in[neighbor[11]]),
            fmedian(in[neighbor[12]], in[neighbor[13]], in[neighbor[14]]),
            fmedian(in[neighbor[15]], in[neighbor[16]], in[neighbor[17]])),
          fmedian(
            fmedian(in[neighbor[18]], in[neighbor[19]], in[neighbor[20]]),
            fmedian(in[neighbor[21]], in[neighbor[22]], in[neighbor[23]]),
            fmedian(in[neighbor[24]], in[neighbor[25]], in[neighbor[26]])));

        // v = max(
        //       max(
        //         max(
        //           max(
        //             max(in[neighbor[0]],in[neighbor[1]]),
        //             max(in[neighbor[2]],in[neighbor[3]])),
        //           max(
        //             max(in[neighbor[4]],in[neighbor[5]]),
        //             max(in[neighbor[6]],in[neighbor[7]]))),
        //         max(
        //           max(
        //             max(in[neighbor[8]],in[neighbor[9]]),
        //             max(in[neighbor[10]],in[neighbor[11]])),
        //           max(
        //             max(in[neighbor[12]],in[neighbor[13]]),
        //             max(in[neighbor[14]],in[neighbor[15]])))),
        //       max(
        //         max(
        //           max(
        //             max(in[neighbor[16]],in[neighbor[17]]),
        //             max(in[neighbor[18]],in[neighbor[19]])),
        //           max(
        //             in[neighbor[20]],
        //             max(in[neighbor[21]],in[neighbor[22]]))),
        //         max(
        //           max(in[neighbor[23]],in[neighbor[24]]),
        //           max(in[neighbor[25]],in[neighbor[26]]))));

        out[xi] = float(v);
      }
    }
  }

  self.curr = itempbuf();
  // printf("max1 -> %p\n", self.buff[self.curr]);
}

void ArFilter::threshold(float min, float max){
  float *data = self.buff[self.curr];
  for(int i=0;i<self.w3;i++){
    if(data[i]<min)data[i]=0;
    if(data[i]>max)data[i]=max;
  }
}
void ArFilter::posterize(int nlevels){
  float *data = self.buff[self.curr];
  for(int i=0;i<self.w3;i++){
    // if(data[i] > 0.1)data[i] = 1;
    // else data[i] = 0;
    // printf("\n%.2f; ", data[i]);
    // data[i] = data[i] * nlevels;
    // printf("-> %.2f; ", data[i]);
    // data[i] = int(data[i]);
    // printf("-> %.2f; ", data[i]);
    // data[i] = data[i] / float(nlevels);
    // printf("-> %.2f; ", data[i]);
    // data[i] = float(int(data[i]*nlevels))/float(nlevels);
  }
}
void ArFilter::normalize(double power){

  NrrdAxisInfo *a = self.a;
  float *data = self.buff[self.curr];
  float *out  = self.buff[itempbuf()];

  float max = 0;
  float min = std::numeric_limits<float>::infinity();
  for(int i=0;i<self.w3;i+=self.w0){
    if(data[i] < min)min=data[i];
    if(data[i] > max)max=data[i];
  }

  max = max-min;
  for(int i=0;i<self.w3;i+=self.w0){
    float r = (data[i]-min)/(max);
    if(power == 1)out[i] = r;
    else out[i] = pow(r,power);
  }
  self.curr = itempbuf();
}
void ArFilter::scale(float s){
  float *data = self.buff[self.curr];
  for(int i=0;i<self.w3;i+=self.w0){
    data[i] *= s;
  }
}
// void ArFilter::median1(){
//   float *in = self.buff[self.curr];
//   float *out = self.buff[itempbuf()];
//   NrrdAxisInfo *a = self.a;

//   int a1 = a[1].size;
//   int a2 = a[2].size;
//   int a3 = a[3].size;

//   float *bufin  = new float[a1*a2*a3];
//   float *bufout = new float[a1*a2*a3];

//   // float *bufin  = (float*) malloc(sizeof(float)*((a->size[1])*(a->size[2])*(a->size[3])));
//   // float *bufout = (float*) malloc(sizeof(float)*((a->size[1])*(a->size[2])*(a->size[3])));
//   int dims[3]  = {a1,a2,a3};
//   int fmin[3]  = {0,0,0};
//   int fsiz[3] = {a1,a2,a3};
//   int fmax[3] = {a1,a2,a3};
//   int l = a1*a2*a3;
//   for(int i=0;i<l;i++){
//     bufin[i] = in[i*2];
//   }
//   median_filter_3D<1>(bufin, dims, bufout, fmin, fsiz, fmax);
//   for(int i=0;i<l;i++){
//     out[i*2] = bufout[i];
//   }
//   free(bufin);
//   free(bufout);
//   self.curr = itempbuf();
// }

void ArFilter::capture(Nrrd *nin){
  if(nin->axis[0].size == self.a0 &&
     nin->axis[1].size == self.a1 &&
     nin->axis[2].size == self.a2){
    memcpy(self.buff[self.curr], nin->data, self.w3 * sizeof(float));
  }
}
void ArFilter::init(Nrrd *nin){
  if(nin == 0){
    self.nrrd = 0;
    self.buff = 0;
    self.nbuf = 0;
    self.curr = 0;
    self.a = 0;

    self.kernel.data = 0;
    self.kernel.temp = 0;
    self.alive = false;
    return;
  }
  if(self.alive){
    destroy();
  }
  self.a = nin->axis;
  
  self.a0=self.a[0].size;
  self.a1=self.a[1].size;
  self.a2=self.a[2].size;
  // self.a2=self.a[3].size;

  // self.w0=1;     // offset to crawl channel
  self.w0=1;               // offset to crawl x
  self.w1=self.a0*self.w0; // offset to crawl y
  self.w2=self.a1*self.w1; // offset to crawl z
  self.w3=self.a2*self.w2; // length of entire dataset.

  self.nbuf = 2;
  self.nrrd = new Nrrd* [self.nbuf];
  self.buff = new float*[self.nbuf];
  self.nrrd[0] = nin;
  self.buff[0] = (float*)nin->data;

  for(int i=1;i<self.nbuf;++i){
    self.nrrd[i] = nrrdNew();
    nrrdCopy(self.nrrd[i],nin);
    self.buff[i] = (float*)(self.nrrd[i]->data);
    // printf("self.buff[i] = %p %p\n", nin->data, self.nrrd[i]->data);
  }
  // printf("filter::init %p %p\n", self.nrrd[0], self.nrrd[1]);
  // printf("filter::init %p %p\n", self.buff[0], self.buff[1]);
  self.curr = 0;
  self.alive = true;
}
Nrrd* ArFilter::commit(Nrrd *nout){
  if(!nout)nout = self.nrrd[0];
  if(0 != nout){
    // printf("memcpy %p %p\n", nout->data, self.nrrd[self.curr]);
    // memset(self.nrrd[self.curr], 0, sizeof(float)*self.w3);
    // printf("commit %p -> %p.\n", self.buff[self.curr], nout);
    memcpy(nout->data, self.buff[self.curr], sizeof(float)*self.w3);
    // printf("done.\n");
    // exit(0);
  }
  return nout;
}
void ArFilter::destroy(){
  // printf("Destroy.\n");
  for(int i=1;i<self.nbuf;++i){
    nrrdNuke(self.nrrd[i]);
  }
  // if(self.kernel.data)delete[] self.kernel.data;
  // if(self.kernel.temp)delete[] self.kernel.temp;
  init(0);
}
int ArFilter::itempbuf(int c){
  c++;
  if(c>=self.nbuf){
    return 0;
  }
  return c;
}
int ArFilter::itempbuf(){
  return itempbuf(self.curr);
}
ArFilter::ArFilter(){
  init(0);
}

ivec3 ArFilter::hill_climb(ivec3 p){
  // printf("%d %d %d; ", p.x, p.y, p.z);
  float *data = self.buff[self.curr];

  if(p.x<1)p.x = 1;
  if(p.y<1)p.y = 1;
  if(p.z<1)p.z = 1;
  if(p.x>self.a0-2)p.x = self.a0-2;
  if(p.y>self.a1-2)p.y = self.a1-2;
  if(p.z>self.a2-2)p.z = self.a2-2;

  float maxv = -1;
  ivec3 maxp;

  int besti = -1;
  int maxi = -1;
  int ii2 = -1;

  bool climbing = true;

  // printf("climbing...");

  while(climbing){
    climbing = false;
    for(int     zz=p.z-1; zz<=p.z+1; zz++){
      for(int   yy=p.y-1; yy<=p.y+1; yy++){
        for(int xx=p.x-1; xx<=p.x+1; xx++){
          // printf("  %d %d %d\n",xx,yy,zz);
          int ii = xx*self.w0 + yy*self.w1 + zz*self.w2;
          float testv = data[ii];
          if(testv > maxv){
            maxv = testv;
            maxp = ivec3(xx,yy,zz);
            // printf("%d %d %d; ", xx, yy, zz);
            climbing = true;
          }
        }
      }
    }
  }
  // printf("\n");
  return maxp;
}
void ArFilter::lapofgaussian_masked(float sigma, int multiply){
  DiscreteKernel kernel = gaussian(sigma, int(sigma*4));
  set_kernel(kernel);
  filter();
  float *orig = new float[self.w3];
  memcpy(orig, self.buff[self.curr], self.w3 * sizeof(float));
  laplacian3d();
  normalize();

  float *out = self.buff[self.curr];
  if(multiply > 0){
    for(int i=0;i<self.w3;i++){
      out[i] *= orig[i];
    }    
  }else if(multiply < 0){
    for(int i=0;i<self.w3;i++){
      if(out[i] > 0) out[i] = orig[i];
    }    
  }
  delete[] orig;
  filter();
  kernel.destroy();

  // highlight(find_maxima());
}
std::vector<glm::ivec3> ArFilter::find_nonzero(){
  std::vector<glm::ivec3> nonzero;
  float *data = self.buff[self.curr];
}
std::vector<glm::ivec3> ArFilter::find_maxima(){
  std::vector<glm::ivec3> maxima;

  float *data = self.buff[self.curr];

  for(int z=0;z<self.a2; z++){
    for(int y=0; y<self.a1; y++){
      for(int x=0; x<self.a0; x++){

        if(x<1 || y<1 || z<1 || x>=self.a0-1 || y>=self.a1-1 || z>=self.a2-1){
          continue;
        }

        int xi = x*self.w0 + y*self.w1 + z*self.w2;
        float p1=0,p2=0,p3=0,
              p4=0,p5=0,p6=0, v;


        v  = data[xi];
        p1 = data[xi-self.w0];
        p2 = data[xi+self.w0];
        p3 = data[xi-self.w1];
        p4 = data[xi+self.w1];
        p5 = data[xi-self.w2];
        p6 = data[xi+self.w2];

        // if(v>0)printf("v: %.4f\n",v);

        if(v>=p1 && v>=p2 && v>=p3 && v>=p4 && v>=p5 && v>=p6 && v != 0 ){
          maxima.push_back(glm::ivec3(x,y,z));
        }
      }
    }
  }
  // printf("found %lu maxima.\n",maxima.size());
  return maxima;
}
void ArFilter::highlight(std::vector<glm::ivec3> points){

  float *buff = self.buff[self.curr];
  for(int i=0;i<self.w3;i+=self.w0){
    buff[i] = buff[i]*0.5f;
  }
  for(int i=0;i<points.size();++i){
    glm::ivec3 v = points[i];
    int xi = v.x*self.w0 + v.y*self.w1 + v.z*self.w2;
    if(xi>=0 && xi<self.w3)buff[xi] = 1.f;
  }
}

// given an input image, find all blobs.
// a blob is a local maximum surrounded by
// its hinterland. ie.:
// a pixel p is in a blob B (centered at pixel b) if 
// p.climb.climb. .... .climb = b.
//
// where p.climb is the adjacent pixel to p with greatest value.
// and b.climb = b.
// we describe this climbing action like a "balloon" which reaches
// its "apex".

typedef unsigned int uint;

struct Balloon{
  uint gradient; // offset to higher balloon (or 0)
  uint volume;   // number of voxels which lead here.
  glm::ivec3 position;
};

struct Point{
  glm::ivec3 p;
  int   i;
};

void ArFilter::show_blobs(int mode){
  float *data = self.buff[self.curr];
  int *labelled = new int[self.w3];
  int nlabels = 0;

  if(!mode)label_blobs(data, labelled);
  else label_blobs_lindeberg(data, labelled);
  // for(int i=0;i<self.w3;i++){
  //   nlabels = max(nlabels, labelled[i]);
  // }
  // printf("nlabels = %d\n", nlabels);
  // float *label_to_value = new float[nlabels];
  // for(int i=0;i<nlabels;i++){
  //   label_to_value[i] = 0;
  //   // for(int j=0;j<self.w3;j++){
  //   //   if(labelled[j] == i){
  //   //     label_to_value[i] = data[j];
  //   //     break;
  //   //   }
  //   // }
  //   label_to_value[i] = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

  // }
  std::map<int, float> labels;
  for(int i=0;i<self.w3;i++){
    int label = labelled[i];
    if(labels.find(label) == labels.end()){
      float v = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      // printf("v=%.8f\n",v);
      labels[label] = v;
      // labels[label] = 0.2f;
    }

    if(labelled[i] >= 0 && data[i]>0.001f){
      // labels[labelled[i]] =0.2f;
      // printf("data[%d] = %d -> %.2f\n", i, labelled[i], labels[labelled[i]]);
      data[i] = labels[labelled[i]];
    }else data[i] = 0;

  }

  printf("nlabels = %d\n", labels.size());
  // delete[] label_to_value;
  delete[] labelled;
}

void ArFilter::label_blobs_lindeberg(float *data, int* labelled){
  long t1 = time(NULL); 
  struct Blob;
  struct Region{
    Region() : pixels(), neighbors(){
      background = false;
      blob = 0;
      id = 0;
      value = 0;
    }
    bool background;              // is this blob background?
    Blob* blob;                   // blob that this region is part of.
    std::unordered_set<int> pixels;         // list of all pixels in this region.
    std::unordered_set<Region*> neighbors;  // list of all neighbors to this region.
    float value;
    int id;
  };
  struct Blob{
    Blob() : regions(){
      cangrow = false;
      id = 0;
    }
    bool cangrow;                    // can this blob grow
    std::vector<Region*> regions;    // list of regions
    int id;                          // unique id label, 0-???
  };

  struct {
    bool operator()(Region *a, Region *b) const{   
        return a->value > b->value;
    }   
  } region_lessthan;

  int *regionmap = new int[self.w3];
  printf("lindeberg.elapsed = %ld %ld %ld\n", (time(NULL)-t1)/3600, (time(NULL)-t1)/60, (time(NULL)-t1));
  std::vector<std::unordered_set<int>> neighbors = label_connected_levelsets(data, regionmap);
  printf("lindeberg.elapsed = %ld %ld %ld\n", (time(NULL)-t1)/3600, (time(NULL)-t1)/60, (time(NULL)-t1));
  // for(std::unordered_set<int> si : neighbors){
  //   printf("%d; ", si.size());
  // }
  // memcpy(labelled, regionmap, sizeof(int)*self.w3);
  // return;
  // printf("\n");

  std::vector<Region*> region_by_value(neighbors.size());
  for(int i=0;i<region_by_value.size();i++)region_by_value[i] = new Region;
  
  std::map<int, Region*> region_by_id;
  
  
  for(int i=0;i<self.w3;i++){
    int label = regionmap[i];
    float value = data[i];
    region_by_value[label]->value = value;
    region_by_value[label]->pixels.insert(i);
    region_by_value[label]->id = label;
    region_by_id[label] = region_by_value[label];
  }
  for(int i=0;i<neighbors.size();i++){
    for(int ni : neighbors[i]){
      region_by_id[i]->neighbors.insert(region_by_id[ni]);
    }
  }
  std::sort(region_by_value.begin(), region_by_value.end(), region_lessthan);

  std::vector<Blob*> blobs;
  printf("lindeberg.elapsed = %ld %ld %ld\n", (time(NULL)-t1)/3600, (time(NULL)-t1)/60, (time(NULL)-t1));
  // printf("regions by value:\n");
  // for(auto x : region_by_value){
  //   printf("%d = %.2f\n", x->id, x->value);
  // }


  // printf("regions by id:\n");
  // for(auto x : region_by_id){
  //   printf("%d = %d, %.2f\n", x.first, x.second->id, x.second->value);
  // }

  for(Region *region : region_by_value){
    if(region->pixels.size() < 1)continue;
    // printf("v = %.3f\n", region->value);
    std::vector<Region*> higher_neighbors;
    bool higher_background_neighbor = false;
    for(Region *n : region->neighbors){
      if(n->value > region->value){
        higher_neighbors.push_back(n);
        if(n->background)higher_background_neighbor = true;
      }
    }
    if(higher_neighbors.size() == 0){
      // case 1: no higher neighbors.
      Blob *blob = new Blob;
      blob->cangrow = true;
      blob->regions.push_back(region);
      blob->id = blobs.size();
      blobs.push_back(blob);

      region->blob = blob;
    }else if(higher_background_neighbor){
      // case 2: higher background neighbor.
      region->background = true;
    }else{
      std::unordered_set<Blob*> neighblobs;
      for(Region *neigh : higher_neighbors){
        if(!neigh->blob){
          fprintf(stderr, "ERROR: REGION WITHOUT BLOB IN CASE 3/4\n");
          exit(0);
        }
        neighblobs.insert(neigh->blob);
      }
      if(neighblobs.size() > 1){
        // case 3: multiple higher blobs. set to background.
        for(Blob *b : neighblobs)b->cangrow = false;
        region->background = true;
      }else{
        // case 4: single higher blob. include region.
        Blob* higher = *neighblobs.begin();
        if(higher->cangrow){
          region->blob = higher;
        }else{
          region->background = true;
        }
      }
    }
  }
  printf("lindeberg.elapsed = %ld %ld %ld\n", (time(NULL)-t1)/3600, (time(NULL)-t1)/60, (time(NULL)-t1));
  // for(int i=0;i<self.w3;i++){        // initialize to -1.
  //   labelled[i] = -1;
  // }
  for(Region *r : region_by_value){
    for(int i : r->pixels){
      labelled[i] = r->blob?r->blob->id:-1;
    }
  }
  // for(int i=0;i<self.w3;i++){        // initialize to -1.
  //   labelled[i] = -1;
  // }
  for(Blob *b : blobs){
    delete b;
  }
  for(Region *r : region_by_value){
    delete r;
  }
  delete[] regionmap;
  printf("}\n");
}

// void ArFilter::label_blobs(float *data, int* labelled){
//   struct regioninfo{
//     regioninfo(){
//       background = false;
//       blob = -1;
//     }
//     bool background;
//     int blob;
//   };
//   struct blobinfo{
//     blob(){
//       cangrow = false;
//     }
//     bool cangrow;
//     std::unordered_set<int> regions;
//   };
//   int nlabels = 0;

//   int *regionmap = new int[self.w3];

//   // neighbors: information about neighbors of regions.
//   // labelled: maps pixel -> region.
//   std::vector<std::unordered_set<int>> neighbors = label_connected_levelsets(data, regionmap);

//   // regions: maps region -> list of pixels
//   std::vector<std::vector<int>> regions(neighbors.size());
//   int nregions = neighbors.size();
//   for(int i=0;i<self.w3;i++){
//     regions[regionmap[i]].push_back(i);
//   }

//   regioninfo* regioninfos = new regioninfo[nregions];
//   for(int i=0;i<nregions;i++){
//     regioninfos[i] = regioninfo();
//   }

//   // list of blobs.
//   std::vector<blob> allblobs;

//   // regionvs: maps region -> value of the particular region levelset.
//   std::vector<std::pair<float, int>> regionvs;
//   std::map<int, float> valueof;
//   for(int r=0;r<regions.size();r++){
//     std::unordered_set<int> pixels = regions[i];
//     if(pixels.size() > 0){
//       int ii = pixels[0];
//       regionvs.push_back(std::make_pair(data[i], r));
//       valueof[r] = data[i];
//     }
//   }

//   sort(regionvs.rbegin(), regionvs.rend());

//   for(int i=0;i<self.w3;i++){        // initialize to -1.
//     labelled[i] = -1;
//     can_grow[i] = false;
//     is_background[i] = false;
//   }

//   int nblobs = 0;

//   for(int i=0;i<nregions;i++){
//     if(regions[i].size() < 1)continue;

//     int   regioni = regionvs[i].second;
//     float regionv = regionvs[i].first;

//     std::unordered_set<int> higher_neighbors;
//     bool higher_background_neighbor = false;
//     for(int n : neighbors[regioni]){
//       float v = valueof[regioni];
//       if(v > regionv){
//         higher_neighbors.insert(i);
//         if(regioninfos[regioni].background){
//           higher_background_neighbor = true;
//         }
//       }
//     }
//     int nhigher_neighbors = higher_neighbors.size();

//     if(nhigher_neighbors == 0){
//       // case 1: no higher neighbor.
//       regioninfos[regioni].blob = nblobs;
//       blobs.push_back(blob());
//       blobs[nblobs].cangrow = true;
//       ++nblobs;
//     }else if(higher_background_neighbor){
//       // case 2: higher neighbor which is background.
//       regioninfos[regioni].background = true;
//     }else{
//       std::unordered_set<int> nblobs;
//       for(int hn : higher_neighbors){
//         nblobs.insert(regioninfos[hn].blob);
//       }
//       if(nblobs.size() > 1){
//         // case 3: higher neighbors, part of different blobs.
//         for(int hn : higher_neighbors){
//           blobs[regioninfos[hn].blob].cangrow = false;
//         }
//         regioninfos[regioni].background = true;
//       }else{
//         // case 4: higher neighbors, part of same blobs.
//         blobs[blobs[0]].cangrow = false;
//       }
//     }
//   }

//   for(int i=)

//   printf("nlabels = %d\n", nlabels);

// }
// void ArFilter::label_blobs(float *data, int* labelled){

//   int nlabels = 0;

//   bool* can_grow = new bool[self.w3];
//   for(int i=0;i<self.w3;i++){        // initialize to -1.
//     labelled[i] = -1;
//     can_grow[i] = false;
//   }


//   std::vector<std::pair<float, int>> pixels;  // [value, index]
//   for(int i=0;i<self.w3;i++){
//     pixels.push_back(std::make_pair(data[i], i));
//   }
//   sort(pixels.rbegin(), pixels.rend()); 

//   for(int i=0;i<self.w3;i++){
//     int index = pixels[i].second;
//     float value = pixels[i].first;

//     int x = (index/self.w0)%self.a0;      // tell the blob where its center is.
//     int y = (index/self.w1)%self.a1;
//     int z = (index/self.w2)%self.a2;

//     // printf("> %.2f\n", x,y,z, value);

//     // float higher_neighborv = 0;
//     int   higher_neighbori = -1;
//     int   higher_neighborl = -1;


//     int higher_neighbors[27];
//     int nhigher_neighbors = 0;
//     bool different_higher_neighbors = false;
//     bool has_higher_neighbor = false;
//     for(int zz=max(z-1,0);zz<=min(z+1,self.a2-1); zz++){
//       for(int yy=max(y-1,0);yy<=min(y+1,self.a1-1); yy++){
//         for(int xx=max(x-1,0);xx<=min(x+1,self.a0-1); xx++){
//           // printf("  %d %d %d\n",xx,yy,zz);
//           int      ii = xx*self.w0 + yy*self.w1 + zz*self.w2;
//           float testv = data[ii];

//           if(testv > value){
//             has_higher_neighbor = true;
//             // higher neighbor...
//             if(labelled[ii] < 0){
//               // case 2: higher neighbor is background.
//               // this pixel must be background.
//               labelled[index] = -1;
//               goto end;
//             }
//             if(higher_neighborl < 0){
//               // first higher neighbor. 
//               higher_neighborl = labelled[ii];  
//             }
//             else{
//               if(higher_neighborl != labelled[ii]){
//                 // case 3: higher neighbors part of different blobs.
//                 // this pixel is background.
//                 different_higher_neighbors = true;
//                 labelled[ii] = -1;
//                 higher_neighbors[nhigher_neighbors] = labelled[ii];
//                 nhigher_neighbors += 1;
//               }
//             }
//           }
//         }
//       }
//     }

//     if(different_higher_neighbors){
//       // case 3: higher neighbors part of different blobs.
//       // stop all higher neighbors from growing.
//       for(int i=0;i<nhigher_neighbors;++i){
//         can_grow[higher_neighbors[i]] = false; 
//       }
//     }
//     else{
//       // case 4: one or more higher neighbors, all part of same blob.
//       if(can_grow[higher_neighborl]){
//         labelled[index] = higher_neighborl;
//       }else{
//         labelled[index] = -1;
//       }
//     }

//     if(!has_higher_neighbor){
//       // case 1: the region has no higher neighbor.
//       labelled[index] = nlabels;
//       can_grow[nlabels] = true;
//       ++nlabels;
//     }

// end:  ;

//     // labelled[index] = x/40;

//   }
//   printf("nlabels = %d\n", nlabels);
//   delete[] can_grow;

// }
void ArFilter::label_blobs(float *data, int* labelled){
  for(int i=0;i<self.w3;i++){        // initialize to -1.
    labelled[i] = -1;
  }

  // breadcumb trail that we leave behind as we look for the max.
  // when we find the max, also update the maxes for each trailing
  // point that we also visited.
  const int max_trail_len = 10000;
  int trail[max_trail_len];
  int trail_length = 0; // length of the trail.

  int x=0,y=0,z=0;

  // hill-climb to determine labels.
  for(int z=0;z<self.a2; z++){
    // if(z%50 == 0)printf("find_blobs progress %d/%d\n",z,self.a2);
    for(int y=0; y<self.a1; y++){
      for(int x=0; x<self.a0; x++){
        int i = x*self.w0 + y*self.w1 + z*self.w2;

        // printf("xy %d %d\n",x,y);
        // starting position is this pixel.
        Point peak;
        peak.p = ivec3(x,y,z);
        peak.i = i;
        // printf("%d %d\n",i, peak.i);

        // printf("(%d %d %d %d)\n",peak.p.x, peak.p.y, peak.p.z, data[peak.i]);
        trail_length = 0;
        for(;;){
          // printf("%d %d\n",i, peak.i);
          // printf("  -> (%d %d %d %d) ",peak.p.x, peak.p.y, peak.p.z, data[peak.i]);
          int   maxi = 0;  // max neighbor index
          float maxv = -1.f;      // max neighbor value
          ivec3 maxp = ivec3(0); // max neighbor coordinates
          // printf("3\n");

          // printf("labelled[%d]\n",peak.i);
          // printf("i=%d %d from %d %d %d\n", i, peak.i, x, y, z);
          if(labelled[peak.i] != -1){
            // we have reached a pixel that already has a label.
            // use this pixel's label as our own (if we keep 
            // climbing we'll reach the same point anyway).
            // printf("labelled[%i]\n",i);
            labelled[i] = labelled[peak.i];
            // printf("4\n");
            break;
          }
          // printf("5\n");

          // add to the trail.
          if(trail_length < max_trail_len){
            trail[trail_length] = peak.i;
            trail_length++;
          }

          // printf("0\n");

          // search for a higher neighbor.
          // do not change the order of traversal.
          // if there are ties, the lowest index is chosen.
          for(int zz=max(peak.p.z-1,0);zz<=min(peak.p.z+1,self.a2-1); zz++){
            for(int yy=max(peak.p.y-1,0);yy<=min(peak.p.y+1,self.a1-1); yy++){
              for(int xx=max(peak.p.x-1,0);xx<=min(peak.p.x+1,self.a0-1); xx++){
                // printf("  %d %d %d\n",xx,yy,zz);
                int      ii = xx*self.w0 + yy*self.w1 + zz*self.w2;
                float testv = data[ii];
                if(testv > maxv){
                  maxv = testv;
                  maxi = ii;
                  maxp = ivec3(xx,yy,zz);
                }
              }
            }
          }
          // printf("1\n");
          if(maxi == peak.i){
            // this is a peak. label it with the index of the maxima.
            // printf("peak at %d = \n",peak.i,peak.i);
            labelled[i] = peak.i;
            break;
          }
          peak.p = maxp;
          peak.i = maxi;
        }
        // printf("2\n");

        // We have successfully labelled this pixel.
        // Now label all the points that we traversed while getting here.
        for(int j=0;j<trail_length;j++){
          labelled[trail[j]] = labelled[i];
        }
      } 
    }
  }
}
std::vector<ScaleBlob*> ArFilter::find_blobs(){
  using namespace glm;

  std::vector<ivec3> maxima;
  // std::vector<glm::ivec3> maxima = find_maxima();
  // printf("maxima: %d\n",maxima.size());

  int *labelled = new int[self.w3];  // mark each point with a label indicated which cell it's part of.

  float *data     = self.buff[self.curr];

  // printf("label blobs.\n");
  label_blobs(data, labelled);


  // now, each pixel knows where its peak is. form a list of unique peaks.

  // float *output = self.buff[self.curr];
  // for(int i=0;i<self.w3;i++){
  //   output[i] = output[i]*3/4;
  // }
  // for(int i=0;i<self.w3;i++){
  //   output[labelled[i]] = 30000;
  // }

  // printf("construct map.\n");
  // form a map, label -> how many voxels are in this blob.
  std::unordered_map<int,int> labels_to_counts;
  for(int i=0;i<self.w3;i+=1){            // form a list of unique labels.
    labels_to_counts[labelled[i]] = 0;
  }
  for(int i=0;i<self.w3;i+=1){            // count how many nonzero pixels have each label.
    if(data[i]>0)
      labels_to_counts[labelled[i]] = labels_to_counts[labelled[i]] + 1;
  }

  // construct a mapping from index to blob.
  // also discard small blobs.
  // also discard negative labels.
  std::unordered_map<int, ScaleBlob*> blobs;
  for ( auto it = labels_to_counts.begin(); it != labels_to_counts.end(); ++it ){
    int index = it->first;
    int count = it->second;
    if(count > 50 && index >= 0){             // impose (arbitrary) minimum blob size, as an optimization.
      blobs[index] = new ScaleBlob();
      float x = (index/self.w0)%self.a0;      // tell the blob where its center is.
      float y = (index/self.w1)%self.a1;
      float z = (index/self.w2)%self.a2;
      blobs[index]->mode = vec3(x,y,z);
      blobs[index]->imode = index;
      blobs[index]->peakvalue = data[index];  // tell the blob the value of its center
    }
  }

  // printf("compute blob statistics.\n");
  // now construct the scaleblobs with data in 3 passes.
  for(int pass=0;pass<2;++pass){
    for(int i=0;i<self.w3;i+=1){
      auto blob = blobs.find(labelled[i]);
      if(blob != blobs.end()){
        float x = (i/self.w0)%self.a0;
        float y = (i/self.w1)%self.a1;
        float z = (i/self.w2)%self.a2;
        if(data[i]>0)blob->second->pass(glm::vec3(x,y,z), (data[i]));
      }
    }
    for (std::pair<int, ScaleBlob*> blob : blobs){
      blob.second->commit();
    }
  }

  // optimize blob centers and shapes using GMM EM algorithm.

  // struct pixel{
  //   float w;    // weight of pixel in blob.
  //   int   i;    // location of pixel.
  //   // data[i] = value of pixel in image.
  // };

  // struct blob_assignment{
  //   // each blob has a list of pixels with weights.
  //   std::vector<pixel> pixels;
  // };

  // std::unordered_map<ScaleBlob*, blob_assignment> assignments;

  // n iterations expectation-maximization algorithm.

  // initialize GMM alpha so each blob as equal weight.
  for (std::pair<int, ScaleBlob*> blob : blobs){
    blob.second->GMM.alpha = 1.f / blobs.size();
  }

  for(int emitr=0;emitr<3;emitr++){
    // expectation step: determine which pixels belong to which blobs.

    // printf("cell positions:\n");
    // for(auto blobitr : blobs){
    //   ScaleBlob *blob = blobitr.second;
    //   printf("%d: %.3f %.3f %.3f\n", blobitr.first, blob->position.x,blob->position.y,blob->position.z);
    // }

    // tree containing all blobs.
    // int blocka0 = self.a0/50 + 2;
    // int blocka1 = self.a1/50 + 2;
    // int blocka2 = self.a2/50 + 2;
    // int blockn  = blocka0 * blocka1 * blocka2;
    // std::vector<ScaleBlob*> *blocks = new std::vector<ScaleBlob*>[blockn];


    // clear weight assignments. add blobs to blocks
    // for(auto blobitr = blobs.begin(); blobitr != blobs.end(); ++blobitr){
    //   blobitr->second->GMM.pixels.clear();
    //   blobitr->second->GMM.weights.clear();
    //   blobitr->second->GMM.nk=0;

    //   int blockx = int(blobitr->second->position.x / 50);
    //   int blocky = int(blobitr->second->position.y / 50);
    //   int blockz = int(blobitr->second->position.z / 50);
    //   int blocki = blockx + blocky*blocka0 + blockz*blocka0*blocka1;
    //   blocks[blocki].push_back(blobitr->second);
    // }

    float N = 0;  // total number of data points.
    printf("d = %d\n", self.w3);


    std::vector<ScaleBlob*> closeblobs;

    float *denominators = new float[self.w3];
    for(int i=0;i<self.w3;i++)denominators[i]=0;
    // memset(denominators, 0, sizeof(float)*self.w3);
    for(int pass=0;pass<2;++pass){
      for(auto blobitr : blobs){
        ScaleBlob *blob = blobitr.second;
        int minx, miny, minz, maxx, maxy, maxz;
        minx = blob->min.x - 20;
        miny = blob->min.y - 20;
        minz = blob->min.z - 20;
        maxx = blob->max.x + 20;
        maxy = blob->max.y + 20;
        maxz = blob->max.z + 20;
        if(minx < 0)minx = 0;
        if(miny < 0)miny = 0;
        if(minz < 0)minz = 0;
        if(maxx >= self.a0-1)maxx = self.a0-1;
        if(maxy >= self.a1-1)maxy = self.a1-1;
        if(maxz >= self.a2-1)maxz = self.a2-1;
        for(int x=minx; x<=maxx; ++x){
          for(int y=miny; y<=maxy; ++y){
            for(int z=minz; z<=maxz; ++z){

              vec3 v(x,y,z);
              int i = x*self.w0 + y*self.w1 + z*self.w2;
              float p = blob->GMM.alpha * blob->gauspdf(v);
              
              if(pass == 0){
                denominators[i] += p;
              }
              
              else{
                float w = (blob->GMM.alpha * blob->gauspdf(v)) / denominators[i];

                if(w > 0.0001){
                  // printf("w=%f\n", w);
                  blob->GMM.pixels.push_back(i);
                  blob->GMM.weights.push_back(w*data[i]);
                  blob->GMM.nk += w * data[i];
                  N += w * data[i];
                }

              }

            }
          }
        }
      }
    }
    delete[] denominators;

    // for(int i=0;i<self.w3;i++){
    //   if(i%100000 == 0) printf("i=%d\n", i);
    //   vec3 v(
    //     (i/self.w0)%self.a0,
    //     (i/self.w1)%self.a1,
    //     (i/self.w2)%self.a2);
    //   // for each pixel
    //   // w(i,k) = numerator / denominator.
    //   float denominator = 0;
    //   // for each blob
    //   // determine ownership of pixel to blob.

    //   closeblobs.clear();

    //   for(int blockx = int(v.x/50)-1; blockx <= int(v.x/50)+1; ++blockx){
    //     for(int blocky = int(v.y/50)-1; blocky <= int(v.y/50)+1; ++blocky){
    //       for(int blockz = int(v.z/50)-1; blockz <= int(v.z/50)+1; ++blockz){
    //         if(blockx<0 || blocky<0 || blockz<0)continue;
    //         int blocki = blockx + blocky*blocka0 + blockz*blocka0*blocka1;
    //         // printf("%d %d %d -> %d / %d\n", blockx, blocky, blockz, blocki, blockn);
    //         closeblobs.insert(closeblobs.end(), blocks[blocki].begin(), blocks[blocki].end());
    //       }
    //     }
    //   }
    //   // printf("l=%d\n", closeblobs.size());
    //   // bspblobs.find_within_distance(closeblobs, v, 2500);

    //   // for(ScaleBlob* blob : closeblobs){
    //   //   // ScaleBlob *blob = blobitr->second;
    //   //   float p = blob->GMM.alpha * blob->gauspdf(v);
    //   //   denominator += p;
    //   // }
    //   // for(ScaleBlob* blob : closeblobs){
    //   //   // ScaleBlob *blob = blobitr->second;
    //   //   float w = (blob->GMM.alpha * blob->gauspdf(v)) / denominator;
    //   //   // printf("w/d = %f/%f\n", w, denominator);
    //   //   if(w > 0.0001){
    //   //     // printf("w=%f\n", w);
    //   //     blob->GMM.pixels.push_back(i);
    //   //     blob->GMM.weights.push_back(w*data[i]);
    //   //     blob->GMM.nk += w * data[i];
    //   //     N += w * data[i];
    //   //   }
    //   // }
    // }
    printf("N = %f\n", N);
    if(N<1e-5)continue;

    // delete[] blocks;

    // maximization step: recalculate mean and covariance.

    for(auto blobitr : blobs){
      ScaleBlob *blob = blobitr.second;
      
      // calculate alpha
      blob->GMM.alpha = blob->GMM.nk / N;
      // printf("alpha = %.3f\n", blob->GMM.alpha);

      // calculate mean, min, max.
      blob->position = dvec3(0);
      float inf = std::numeric_limits<float>::infinity();
      blob->min = vec3(inf, inf, inf);
      blob->max = vec3(0);
      for(int i=0;i<blob->GMM.pixels.size();i++){
        vec3 v(
          (blob->GMM.pixels[i]/self.w0)%self.a0,
          (blob->GMM.pixels[i]/self.w1)%self.a1,
          (blob->GMM.pixels[i]/self.w2)%self.a2);
        blob->position += v * blob->GMM.weights[i];
        if(v.x<blob->min.x)blob->min.x = v.x;
        if(v.y<blob->min.y)blob->min.y = v.y;
        if(v.z<blob->min.z)blob->min.z = v.z;
        if(v.x>blob->max.x)blob->max.x = v.x;
        if(v.y>blob->max.y)blob->max.y = v.y;
        if(v.z>blob->max.z)blob->max.z = v.z; 
      }
      blob->position /= blob->GMM.nk;

      // calculate covariance
      blob->shape = mat3x3(0);
      for(int i=0;i<blob->GMM.pixels.size();i++){
        vec3 v(
          (blob->GMM.pixels[i]/self.w0)%self.a0,
          (blob->GMM.pixels[i]/self.w1)%self.a1,
          (blob->GMM.pixels[i]/self.w2)%self.a2);
        v -= blob->position;
        blob->shape[0][0] += v.x*v.x*blob->GMM.weights[i];
        blob->shape[0][1] += v.x*v.y*blob->GMM.weights[i];
        blob->shape[0][2] += v.x*v.z*blob->GMM.weights[i];
        blob->shape[1][1] += v.y*v.y*blob->GMM.weights[i];
        blob->shape[1][2] += v.y*v.z*blob->GMM.weights[i];
        blob->shape[2][2] += v.z*v.z*blob->GMM.weights[i];
      }
      blob->shape /= blob->GMM.nk;

      blob->npass = 1;
      blob->commit();

    }



  }
  // printf("list blobs:\n");
  // for (std::pair<int, ScaleBlob*> blob : blobs){
  //   blob.second->print();
  // }

  // std::vector<ivec3>       positions;
  std::vector<ScaleBlob*> output;

  for (std::pair<int, ScaleBlob*> blob : blobs){
    // another optimization: only add blobs with volume > 4.
    if(blob.second->detCov > 10.f && blob.second->n > 30 ){
      output.push_back(blob.second);
    }
    // positions.push_back(ivec3(blob.second->position));
  }

  // construct a sorted list of unique labels from this set.
  // std::vector<int> labels;
  // labels.assign( set.begin(), set.end() );
  // sort( labels.begin(), labels.end() );

  // construct a list of counts.
  // std::vector<int> counts(labels.size());

  // std::vector<ivec3>
  // normalize(0.1);
  // highlight(positions);
  // printf("%lu blobs.\n",output.size());
  delete[] labelled;
  // printf("done.\n");
  return output;
}

void ArFilter::print(){

  int channel = 0;  
  int buckets[10];
  for(int i=0;i<10;++i)buckets[i] = 0;

  printf("axes %d %d %d %d\n",self.a0,self.a0,self.a1,self.a2);

  NrrdAxisInfo *a = self.a;
  float *data = self.buff[self.curr];

  double max = 0;
  double min = std::numeric_limits<double>::infinity();
  for(int i=channel;i<self.w3;i+=2){
    if(data[i] < min)min=data[i];
    if(data[i] > max)max=data[i];
  }
  max = max-min;
  for(int i=channel;i<self.w3;i+=2){
    double r = double(data[i]-min)/double(max);
    int bucket = (int)(r*10.0);
    if(bucket<0)bucket=0;
    if(bucket>9)bucket=9;
    buckets[bucket]++;
  }
  printf("histogram:\n ");
  for(int i=0;i<10;i++){
    printf("%4d, ",buckets[i]);
  }
  printf("\n");
  printf("minmax: %.1f %.1f\n",min,min+max);
}

void ArFilter::clear(){
  float *data = self.buff[self.curr];
  for(int i=0;i<self.w3;i++)data[i] = 0;
}
void ArFilter::rasterlineadd(vec3 a, vec3 b, float va, float vb){
  float *data = self.buff[self.curr];
  int i;

  vec3 v = b-a;
  float len = length(b-a);
  vec3 step = v/len;

  for(int j=0;j<len;++j){
    a += step;
    i = int(a.x)*self.w0 + int(a.y)*self.w1 + int(a.z)*self.w2;
    // printf(". %.1f %.1f %.1f\n", a.x, a.y, a.z);
    data[i] = va + (float(j)/len)*(vb-va);
  }
}
void ArFilter::color_blobs(std::vector<ScaleBlob*> blobs, float color){
  using namespace glm;
  float *data = self.buff[self.curr];

  // iterate through blobs.
  for(auto blob = blobs.begin(); blob != blobs.end(); ++blob){
    ScaleBlob *sb = *blob;
    float minx = sb->min.x;
    float miny = sb->min.y;
    float minz = sb->min.z;
    float maxx = sb->max.x;
    float maxy = sb->max.y;
    float maxz = sb->max.z;

    // iterate through pixels for each blob.
    for(float x=minx; x<=maxx; ++x){
      for(float y=miny; y<=maxy; ++y){
        for(float z=minz; z<=maxz; ++z){
          int i = int(x)*self.w0 + int(y)*self.w1 + int(z)*self.w2;
          float v = sb->cellpdf(vec3(x,y,z));
          if(std::isfinite(v) && v > 0.0001f){
            float orig = data[i] - 2.f*int(data[i]/2.f);
            data[i] = color + max(orig,v);
          }
        }
      }
    }
  }
}
struct thread_drawblobs_info{
  ArFilter *filter;
  std::vector<ScaleBlob*> *blobs;
  float *data;
  float *lock;
  int blobmin;          // render window min value.
  int blobmax;          // render window max value.
  const char *mode;     // either 'q'=quick or 'g'=gaussian
                        // mode[2] either '+' or 'm'
};

static void* t_draw_blobs(void* vinfo){
  thread_drawblobs_info *info = (thread_drawblobs_info*)vinfo;
  float *data = info->data;
  float *lock = info->lock;
  int blobmin = info->blobmin;
  int blobmax = info->blobmax;
  const char *mode   = info->mode;
  char mode0 = mode[0];
  ArFilter *filter = info->filter;
  std::vector<ScaleBlob*> &blobs = *info->blobs;

  // printf("t_draw_blobs %d - %d\n", blobmin, blobmax);
  for(int bi=blobmin; bi<blobmax; ++bi){
    ScaleBlob *sb = blobs[bi];
            // printf("dot");
    if(sb->n < 3)continue;
    int minx, miny, minz, maxx, maxy, maxz;
    if(sb->model.type != '_'){
      minx = sb->model.min.x;
      miny = sb->model.min.y;
      minz = sb->model.min.z;
      maxx = sb->model.max.x;
      maxy = sb->model.max.y;
      maxz = sb->model.max.z;
      if(  sb->model.type == 'g'
        || sb->model.type == 'f')mode0 = sb->model.type;
    }else{
      minx = sb->min.x;
      miny = sb->min.y;
      minz = sb->min.z;
      maxx = sb->max.x;
      maxy = sb->max.y;
      maxz = sb->max.z;
    }
    if(minx < 0)minx = 0;
    if(miny < 0)miny = 0;
    if(minz < 0)minz = 0;
    if(maxx >= filter->self.a0)maxx = filter->self.a0;
    if(maxy >= filter->self.a1)maxy = filter->self.a1;
    if(maxz >= filter->self.a2)maxz = filter->self.a2;
    float v;
    // printf("minmax %d %d %d %d %d %d\n",minx, maxx, miny, maxy, minz, maxz);
    for(int x=minx; x<=maxx; ++x){
      for(int y=miny; y<=maxy; ++y){
        for(int z=minz; z<=maxz; ++z){
          int i = (x)*filter->self.w0 + (y)*filter->self.w1 + (z)*filter->self.w2;
          if(i<0)continue;
          if(i>=filter->self.w3)continue;
          // float 
          if(mode0 == 'g'){
            v = sb->pdf(vec3(x,y,z));
          }
          else if(mode0=='q'){
            v = sb->cellpdf(vec3(x,y,z));
          }
          else if(mode0=='f'){
            v = sb->cellerf(vec3(x,y,z));
          }
          else if(mode0=='.'){
            v = sb->celldot(vec3(x,y,z));
          }
          else if(mode0=='l'){
            v = sb->ellipsepdf(vec3(x,y,z));
          }
          else if(mode0=='o'){
            v = sb->outlinepdf(vec3(x,y,z));
          }
          else if(mode0=='m'){
            v = sb->generalized_multivariate_gaussian_pdf(vec3(x,y,z));
          }
          else if(mode0=='M'){
            v = sb->gmmpdf(vec3(x,y,z));
          }
          // float v = 0.5f;
          // if(std::isfinite(v)){
          if(mode[1] == '+'){
            data[i] = (data[i] + v);
            if(data[i]>1)data[i] = 1;
          }
          else data[i] = max(data[i],v);
            // printf("(%.2f %.2f %.2f) -> %.f", x,y,z,v);
          // }
        }
      }
    }
  }
  return 0;
}

void ArFilter::draw_blobs(std::vector<ScaleBlob*> blobs, const char *mode){
  // printf("draw blobs.");
  using namespace glm;

  const int nthreads = 12;
  float *lock  = self.buff[itempbuf()];
  float *data  = self.buff[self.curr];

  pthread_t threads[nthreads];
  thread_drawblobs_info info[nthreads];


  for(int i=0;i<nthreads;i++){
    info[i].blobs = &blobs;
    info[i].lock = lock;
    info[i].data = data;
    info[i].filter = this;

    info[i].blobmin = (i*blobs.size())/nthreads;
    info[i].blobmax = ((i+1)*(blobs.size()))/nthreads;
    info[i].mode    = mode;
  }

  for(int i=0;i<nthreads;++i){
    if(pthread_create(threads+i, NULL, t_draw_blobs, info+i)){
      fprintf(stderr, "error creating render thread.\n");
      exit(0);
    }
  }
  for(int i=0;i<nthreads;++i){
    if(int err = pthread_join(threads[i], 0)){
      fprintf(stderr, "error joining render thread. %d\n", err);
      exit(0);
    }
  }


  // tick("draw blobs");
  // printf("done.");

  /////////////////////////////////////////////////////////////////////

  // ** Normalizing shouldn't be necessary, because cellpdf < 1 ** .
  // // get max value.
  // float max = 0;
  // for(int i=0;i<self.w3;i+=2){
  //   if(data[i] > max)max=data[i];
  // }

  // // divide by max value.
  // for(int i=0;i<self.w3;i+=2){
  //   data[i] = (data[i])/(max);
  // }

  ////////////////////////////////////////////////////////////////////

  // scale(0.2f);
  // for(auto blob = blobs.begin(); blob != blobs.end(); ++blob){
  //   ScaleBlob *sb = *blob;
  //   if(sb->parent){
  //     rasterlineadd(vec3(sb->position), vec3(sb->parent->position), 1.1f, 1.5f);
  //   }
  //   // for(ScaleBlob *succ : sb->succ){
  //   //   rasterlineadd(vec3(sb->position), vec3(succ->position), 1.1f, 1.5f);
  //   // }
  // }
  // for(int z=0;z<self.a2; z++){
  //   if(z%50 == 0)printf("draw_blobs progress %d/%d\n",z,self.a2);
  //   for(int y=0; y<self.a1; y++){
  //     for(int x=0; x<self.a0; x++){
  //       int i = x*self.w0 + y*self.w1 + z*self.w2;
  //       vec3 p(x,y,z);
  //       float v= 0.f;
  //       for ( auto blob = blobs.begin(); blob != blobs.end(); ++blob ){
  //         float vv = (*blob)->pdf(p);
  //         if(std::isfinite(vv)){
  //           v += vv;
  //         }
  //       }
  //       data[i] = v;
  //     }
  //   }
  // }
  // normalize(1.f);
  return;

  // if(highlight){
  //   normalize(0.25);
  //   // highlight centers.
  //   float *buff = self.buff[self.curr];
  //   for(int i=0;i<self.w3;i+=self.w0){
  //     buff[i] = buff[i]*0.95f;
  //   }
  //   for(int i=0;i<blobs.size();++i){
  //     glm::ivec3 v(blobs[i]->position);
  //     int xi = v.x*self.w0 + v.y*self.w1 + v.z*self.w2;
  //     if(xi>=0 && xi<self.w3)buff[xi] = 1.f;
  //   }
  // }
}

static inline void fppswap(float **a, float **b){
  float *temp = *a;
  *b          = *a;
  *a          = temp;
}
static inline float sq(float x){
  return x*x;
}
static int median(std::vector<int> &v){
  size_t n = v.size() / 2;
  nth_element(v.begin(), v.begin()+n, v.end());
  return v[n];
}
static int mean(std::vector<int> &v){
  int sum = 0;
  for(int x : v){
    sum += x;
  }
  return sum/v.size();
}
int ArFilter::count_blobs(){
  float *data   = self.buff[self.curr];
  int *labelled = new int[self.w3];
  label_blobs(data, labelled);


  // now, each pixel knows where its peak is. form a list of unique peaks.

  // float *output = self.buff[self.curr];
  // for(int i=0;i<self.w3;i++){
  //   output[i] = output[i]*3/4;
  // }
  // for(int i=0;i<self.w3;i++){
  //   output[labelled[i]] = 30000;
  // }

  // form a map, label -> how many voxels are in this blob.
  std::map<int,int> labels_to_counts;
  for(int i=0;i<self.w3;i+=1){            // form a list of unique labels.
    labels_to_counts[labelled[i]] = 0;
  }
  for(int i=0;i<self.w3;i+=1){            // count how many pixels have each label.
    labels_to_counts[labelled[i]] = labels_to_counts[labelled[i]] + 1;
  }
  std::vector<int> counts;
  for(auto x : labels_to_counts){
    if(x.second > 7000) counts.push_back(x.second);
  }
  delete[] labelled;
  return counts.size();
}

std::vector<std::unordered_set<int>> ArFilter::label_connected_levelsets(float *data, int* labels){
  long t1 = time(NULL);
  // float *data   = self.buff[self.curr];
  std::vector<std::unordered_set<int>> neighbors;
  for(int i=0;i<self.w3;i++){
    labels[i] = 0;
  }
  int i;
  int label = 0;
  float value = -1;
  std::queue<int> traverse;   // queue to mark entire connected component.
  printf("levelsets.elapsed = %ld %ld %ld\n", (time(NULL)-t1)/3600, (time(NULL)-t1)/60, (time(NULL)-t1));
  for(int i=0;i<self.w3;i++){
    if(labels[i] != 0)continue; // already labeled.
    if(data[i] == 0)continue;   // part of background. ignore.
    traverse = std::queue<int>();
    traverse.push(i);
    // find entire connected component.
    ++label;           // label for levelset region.
    value = data[i];   // value of levelset region.
    int ct = 0;
    while(!traverse.empty()){
      // printf(" size = %d\n", traverse.size());
      int curr = traverse.front();
      // printf("c[%d] = %d\n", curr, labels[curr]);
      traverse.pop();
      // printf("pop!\n");

      // where are we?

      bool valid = (curr>=0 && curr<self.w3);

      // bool valid = x>=0 && y>=0 && z>=0 && 
      //       x<self.a0 && y<self.a1 && z<self.a2 && 
      //       curr>=0 && curr<self.w3;

      // if(curr >= self.w3){
      //   fprintf(stderr,"bad! %d %d %d %d\n",curr, x,y,z);
      //   exit(0);
      // }

      // printf("%d %d %d %d\n", curr, int(valid), int(labels[curr]), int(data[curr]!=0));

      if( valid && !labels[curr] && data[curr]!=0 && fabs(data[curr] - value) < 0.02f){
        int x = (curr/self.w0)%self.a0;
        int y = (curr/self.w1)%self.a1;
        int z = (curr/self.w2)%self.a2;
        if(x<1 || y<1 || z<1 || x>=self.a0-1 || x>=self.a1-1 || x>=self.a2-1)continue;
        // unlabeled and a real point and nonzero and correct value
        // printf("labels[%d] = %d\n", curr, label);
        labels[curr] = label;

        traverse.push(curr + -1*self.w0 + -1*self.w1 + -1*self.w2);
        traverse.push(curr + -1*self.w0 + -1*self.w1 +  0*self.w2);
        traverse.push(curr + -1*self.w0 + -1*self.w1 +  1*self.w2);
        traverse.push(curr + -1*self.w0 +  0*self.w1 + -1*self.w2);
        traverse.push(curr + -1*self.w0 +  0*self.w1 +  0*self.w2);
        traverse.push(curr + -1*self.w0 +  0*self.w1 +  1*self.w2);
        traverse.push(curr + -1*self.w0 +  1*self.w1 + -1*self.w2);
        traverse.push(curr + -1*self.w0 +  1*self.w1 +  0*self.w2);
        traverse.push(curr + -1*self.w0 +  1*self.w1 +  1*self.w2);

        traverse.push(curr +  0*self.w0 + -1*self.w1 + -1*self.w2);
        traverse.push(curr +  0*self.w0 + -1*self.w1 +  0*self.w2);
        traverse.push(curr +  0*self.w0 + -1*self.w1 +  1*self.w2);
        traverse.push(curr +  0*self.w0 +  0*self.w1 + -1*self.w2);
     // traverse.push(curr +  0*self.w0 +  0*self.w1 +  0*self.w2);
        traverse.push(curr +  0*self.w0 +  0*self.w1 +  1*self.w2);
        traverse.push(curr +  0*self.w0 +  1*self.w1 + -1*self.w2);
        traverse.push(curr +  0*self.w0 +  1*self.w1 +  0*self.w2);
        traverse.push(curr +  0*self.w0 +  1*self.w1 +  1*self.w2);

        traverse.push(curr +  1*self.w0 + -1*self.w1 + -1*self.w2);
        traverse.push(curr +  1*self.w0 + -1*self.w1 +  0*self.w2);
        traverse.push(curr +  1*self.w0 + -1*self.w1 +  1*self.w2);
        traverse.push(curr +  1*self.w0 +  0*self.w1 + -1*self.w2);
        traverse.push(curr +  1*self.w0 +  0*self.w1 +  0*self.w2);
        traverse.push(curr +  1*self.w0 +  0*self.w1 +  1*self.w2);
        traverse.push(curr +  1*self.w0 +  1*self.w1 + -1*self.w2);
        traverse.push(curr +  1*self.w0 +  1*self.w1 +  0*self.w2);
        traverse.push(curr +  1*self.w0 +  1*self.w1 +  1*self.w2);
        
        // printf("push6\n");
        // traverse.push(curr+self.w0);
        // traverse.push(curr-self.w0);
        // traverse.push(curr+self.w1);
        // traverse.push(curr-self.w1);
        // traverse.push(curr+self.w2);
        // traverse.push(curr-self.w2);

        // data[curr] = r;
        ++ct;

      }
    }
    // break;
    // printf("found cpt of size %d = %d\n", label, ct);
  }
  printf("levelsets.elapsed = %ld %ld %ld\n", (time(NULL)-t1)/3600, (time(NULL)-t1)/60, (time(NULL)-t1));

  while(neighbors.size() <= label)neighbors.push_back(std::unordered_set<int>());
  for(int i=0;i<self.w3;i++){
    int x = (i/self.w0)%self.a0;
    int y = (i/self.w1)%self.a1;
    int z = (i/self.w2)%self.a2;
    if(x<1 || y<1 || z<1 || x>=self.a0-1 || y>=self.a1-1 || z>=self.a2-1)continue;
    neighbors[labels[i]].insert(labels[i-1]);

    neighbors[labels[i]].insert(labels[i + -1*self.w0 + -1*self.w1 + -1*self.w2]);
    neighbors[labels[i]].insert(labels[i + -1*self.w0 + -1*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i + -1*self.w0 + -1*self.w1 +  1*self.w2]);
    neighbors[labels[i]].insert(labels[i + -1*self.w0 +  0*self.w1 + -1*self.w2]);
    neighbors[labels[i]].insert(labels[i + -1*self.w0 +  0*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i + -1*self.w0 +  0*self.w1 +  1*self.w2]);
    neighbors[labels[i]].insert(labels[i + -1*self.w0 +  1*self.w1 + -1*self.w2]);
    neighbors[labels[i]].insert(labels[i + -1*self.w0 +  1*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i + -1*self.w0 +  1*self.w1 +  1*self.w2]);

    neighbors[labels[i]].insert(labels[i +  0*self.w0 + -1*self.w1 + -1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  0*self.w0 + -1*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i +  0*self.w0 + -1*self.w1 +  1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  0*self.w0 +  0*self.w1 + -1*self.w2]);
 // neighbors[labels[i]].insert(labels[i +  0*self.w0 +  0*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i +  0*self.w0 +  0*self.w1 +  1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  0*self.w0 +  1*self.w1 + -1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  0*self.w0 +  1*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i +  0*self.w0 +  1*self.w1 +  1*self.w2]);

    neighbors[labels[i]].insert(labels[i +  1*self.w0 + -1*self.w1 + -1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  1*self.w0 + -1*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i +  1*self.w0 + -1*self.w1 +  1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  1*self.w0 +  0*self.w1 + -1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  1*self.w0 +  0*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i +  1*self.w0 +  0*self.w1 +  1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  1*self.w0 +  1*self.w1 + -1*self.w2]);
    neighbors[labels[i]].insert(labels[i +  1*self.w0 +  1*self.w1 +  0*self.w2]);
    neighbors[labels[i]].insert(labels[i +  1*self.w0 +  1*self.w1 +  1*self.w2]);

    // for(int zz=max(z-1,0);zz<=min(z+1,self.a2-1); zz++){
    //   for(int yy=max(y-1,0);yy<=min(y+1,self.a1-1); yy++){
    //     for(int xx=max(x-1,0);xx<=min(x+1,self.a0-1); xx++){
    //       int ii = xx*self.w0 + yy*self.w1 + zz*self.w2;
    //       if(i == ii)continue;
    //       neighbors[labels[i]].insert(labels[ii]);
    //     }
    //   }
    // }
  }

  return neighbors;
  // delete[] labels;
  // return label;
}

int ArFilter::count_connected_components(){
  float *data   = self.buff[self.curr];
  int *labels = new int[self.w3];
  for(int i=0;i<self.w3;i++){
    labels[i] = 0;
  }
  int i;
  int label = 0;
  for(int i=0;i<self.w3;i++){
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    // r *= 3;
    // i = x*self.w0 + y*self.w1 + z*self.w2;
    if(labels[i] != 0)continue; // already labeled.
    if(data[i] == 0)continue;   // part of background. ignore.
    std::queue<int> traverse;   // queue to mark entire connected component.
    traverse.push(i);
    
    // find entire connected component.
    ++label;  // label for connected component.
    int ct = 0;
    while(!traverse.empty()){
      int curr = traverse.front();
      traverse.pop();

      // where are we?
      int x = (curr/self.w0)%self.a0;
      int y = (curr/self.w1)%self.a1;
      int z = (curr/self.w2)%self.a2;

      bool valid = x>=0 && y>=0 && z>=0 && 
            x<self.a0 && y<self.a1 && z<self.a2;

      // printf("%d %d %d %d\n", curr, int(valid), int(labels[curr]), int(data[curr]!=0));

      if( valid && !labels[curr] && data[curr]!=0 ){
        // unlabeled and a real point and nonzero
        // printf("labels[%d] = %d\n", curr, label);
        labels[curr] = label;
        
        traverse.push(curr+self.w0);
        traverse.push(curr-self.w0);
        traverse.push(curr+self.w1);
        traverse.push(curr-self.w1);
        traverse.push(curr+self.w2);
        traverse.push(curr-self.w2);

        data[curr] = r;
        ++ct;

      }
    }
    // break;
    // printf("found cpt of size %d = %d\n", label, ct);
  }

  delete[] labels;
  return label;
}


ScaleBlob* ArFilter::compute_blob_tree(){
  printf("compute_blob_tree...");
  BSPTree<ScaleBlob> bspblobs(0,10,vec3(-1.f,-1.f,-1.f),vec3(self.a0+1.f,self.a1+1.f,self.a2+1.f)); // safe against any rounding errors.

  DiscreteKernel kernel; 
  float sigma = .7f;
  float scale = 0.f;

  float *const iscalec = new float[self.w3];   // data for scaled image.
  float *iscale = iscalec;                     // so we can manipulate this pointer.

  std::vector<ScaleBlob*> blobs;
  float *temp;

  int index_of_laplacian;

  // printf("threshold 0.1f;\n");
  // threshold(0.1f, 1.f);
  max1();
  // we want to compute the gaussian of the image at various scales.
  // use the fact that the composition of gaussians is a gaussian with
  //  c^2 = a^2 + b^2.
  while(scale < 6.f){

    // printf("filter stack: %p %p. curr=%d.\n", self.buff[0], self.buff[1], self.curr);
    scale = sqrt(scale*scale + sigma*sigma);  // compute scale factor with repeated gaussians.
    // printf("gaussian %.2f -> scale %.2f",sigma, scale);
    kernel = gaussian(sigma, int(sigma*4));   // create gaussian kernel with new sigma.
    set_kernel(kernel);                       //
    filter();                                 // compute gaussian blur, store in self.buff[curr].
    printf("scale = %.2f; ", scale);
    float **blurred = self.buff + self.curr;  // address of blurred image.
    // printf("  laplacian.. ");
    // pass();
    laplacian3d(1);                           // compute laplacian, store in self.buff[curr+1].
    // threshold(0.01,1);
    // fppswap(blurred, &iscale);                 // remove the blurred image from the swap chain
                                              // and store in iscale, so that it's not overwritten
                                              // by successive filter operations.
    for(int i=0;i<self.w3;i++){
      if(self.buff[self.curr][i] > 0){
        self.buff[self.curr][i] = (*blurred)[i];
      }
    }
    // index_of_laplacian = self.curr;

    temp     = iscale;
    iscale   = *blurred;
    *blurred = temp;

    // printf("swap.\n");
    // printf("filter stack: %p %p. curr=%d.\n", self.buff[0], self.buff[1], self.curr);
    // printf("  normalize.. ");

    // float *curr = self.buff[self.curr];
    // for(int i=0;i<self.w3;i++){
    //   if(curr[i] > 0){
    //     curr[i] = iscale[i];
    //   }
    // }
    // filter();
    normalize();
    // printf(". find blobs.\n");
    std::vector<ScaleBlob*> bigblobs = find_blobs();
    printf("found %d blobs. ", bigblobs.size());
    for(ScaleBlob *sb : bigblobs){
      sb->scale = scale;
      sb->peakvalue = temp[sb->imode];
      bspblobs.insert(sb, sb->mode);
    }
    // printf("connect..\n");
    for(ScaleBlob *x : blobs){
      // ScaleBlob *closest = 0;
      // float      min_dist = -1;
      // for(ScaleBlob *y : bigblobs){
      //   float dist = sq(x->mode.x - y->mode.x) + sq(x->mode.y - y->mode.y) + sq(x->mode.z - y->mode.z);
      //   if(!closest || dist < min_dist){
      //     min_dist = dist;
      //     closest = y;
      //   }
      // }

      if(!x)continue;

      ScaleBlob *closest  = bspblobs.find_closest(x->mode);
      x->parent = closest;
      closest->children.push_back(x);
    }
    bspblobs.clear();

    blobs = bigblobs;

    // prune the tree so that blobs never only have one child.
    // for(int i=0;i<blobs.size();++i){
    //   if(blobs[i]->children.size() == 1){
    //     ScaleBlob *only_child = blobs[i]->children[0];
    //     blobs[i] = only_child;
    //     only_child->parent = 0;
    //     // delete *blobs[i];
    //   }
    // }

    sigma += 0.2f;
    // sigma += 1;
    // sigma = sqrt(2*scale + 0.5);
    kernel.destroy();

    temp = self.buff[self.curr];
    self.buff[self.curr] = iscale;
    iscale = temp;

    // fppswap(&iscale, self.buff+self.curr);    // put the gaussian blurred image back at self.buff[curr]
                                              // for use in the next filter operation.
    // printf(".done.\n");
  }
  printf("..done.\n");

  // undo our pointer shenanigans.
  // remove the buffer that we just created
  // from the swap chain.
  for(int i=0;i<self.nbuf;++i){
    if(self.buff[i] == iscalec){
      self.buff[i] = iscale;
    }
  }
  // printf("filter stack: %p %p. iscalec=%p.\n", self.buff[0], self.buff[1], iscalec);
  --self.curr;
  if(self.curr<0)self.curr = self.nbuf-1;
  normalize();
  // self.curr = index_of_laplacian;
  // laplacian3d();
  delete[] iscalec;
  ScaleBlob *frame = new ScaleBlob();
  frame->children = blobs;
  frame->parent = 0;
  return frame;
}

BSPTree<ScaleBlob> ArFilter::get_bsp(int depth){
  return BSPTree<ScaleBlob>(0, depth, vec3(-1.f,-1.f,-1.f),vec3(self.a0+1.f,self.a1+1.f,self.a2+1.f));
}

void ArFilter::difference_image(Nrrd* x){
  float *data1 = self.buff[self.curr];
  float *data2 = (float*) x->data;
  for(int i=0;i<self.w3;i+=self.w0){
    data1[i] = fabs(data1[i] - data2[i]);
  }
}

#undef pi
#undef gaus
#undef gausd2