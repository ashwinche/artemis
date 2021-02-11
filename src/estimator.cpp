#include "estimator.h"
#include <limits>
#include <random>
#include <set>
#include <functional>


static float randf(){
  return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
}
static float randf(float min, float max){
  return min + randf()*(max-min);
}
ScaleBlob Estimator::fit(Nrrd* source, ScaleBlob* in){
  ScaleBlob blob;
  blob.position = in->position;
  blob.min = in->min;
  blob.max = in->max;
  blob.n = in->n;
  blob.shape = in->shape;
  blob.mode = in->mode;
  blob.invCov = in->invCov;

  int a0 = source->axis[0].size;
  int a1 = source->axis[1].size;
  int a2 = source->axis[2].size;


  float bestalpha = 0;
  float bestbeta = 0;
  float bestkappa = 0;
  float besterror = std::numeric_limits<float>::infinity();

  float error = 0;
  std::vector<ivec3> indices;
  int n = 0;
  printf("indices..");
  int xmin = blob.position.x - (blob.position.x - blob.min.x)*2;
  int ymin = blob.position.y - (blob.position.y - blob.min.y)*2;
  int zmin = blob.position.z - (blob.position.z - blob.min.z)*2;

  int xmax = blob.position.x + (blob.max.x - blob.position.x)*2;
  int ymax = blob.position.y + (blob.max.y - blob.position.y)*2;
  int zmax = blob.position.z + (blob.max.z - blob.position.z)*2;

  if(xmin<0)xmin =0;
  if(ymin<0)ymin =0;
  if(zmin<0)zmin =0;

  if(xmax>=a0)xmax = a0 - 1;
  if(ymax>=a1)ymax = a1 - 1;
  if(zmax>=a1)zmax = a1 - 1;

  for (int xx=xmin; xx < xmax; xx++){
    for(int yy=ymin; yy < ymax; yy++){
      for(int zz=zmin; zz < zmax; zz++){
         vec3 p(xx,yy,zz);
         if(blob.ellipsepdf(p) > 0){
          indices.push_back(ivec3(xx,yy,zz));
          n = n+1;
        }
      }
    }
  }
  if(n == a0*a1*a2){
    indices = std::vector<ivec3>();
    indices.push_back(ivec3(blob.position.x, blob.position.y, blob.position.z));
  }
  // std::function<void(vec3)> pdf(std::bind(&ScaleBlob::generalized_multivariate_gaussian_pdf, blob, _1));
  // printf("max = %d\n",a0*a1*a2);

  // char model_type = 'g';
  // auto pdf = [&blob](vec3 v) { return blob.generalized_multivariate_gaussian_pdf(v); };

  char model_type = 'f';
  auto pdf = [&blob](vec3 v) { return blob.cellerf(v); };

  printf("estimate.. n=%d, ", n);
  float beta = 1.f;
  for(int i=0;i<200;i++){
    // printf("%d, ", i);
    float beta  = randf(1.f,3.f);
    float kappa = randf(0.1,1.f);
    float alpha = randf(0.2, 0.5);

  // for(float beta = 1.f; beta < 6.f; beta *= 1.2f){
  //   // printf("beta = %.2f. ", beta);
  //   for(float kappa = 0.1; kappa < 0.95; kappa += 0.05f){
  //     for(float alpha=0.08f; alpha < 0.5f; alpha*= 1.1){
        // float r1 = 
        // float beta = 1.f;
        // for each beta, kappa, alpha:
        blob.model.kappa = kappa;
        blob.model.alpha = alpha;
        blob.model.beta  = beta;
        error = 0;
        // printf("hi . ");
        for(ivec3 v : indices){
          vec3 p(v);
          int xx = v.x;
          int yy = v.y;
          int zz = v.z;
          int i = xx + source->axis[0].size*yy + source->axis[0].size*source->axis[1].size*zz;
          float diff = ((float*)source->data)[i] - pdf(p);
          error += diff*diff;
        }
        // error /= n;
        // printf("  error = %.5f\n", error);
        if(error<besterror){
          besterror = error;
          bestalpha = alpha;
          bestbeta  = beta;
          bestkappa = kappa;
        }
    //   }
    // }
  }
  blob.model.alpha = bestalpha;
  blob.model.beta = bestbeta;
  blob.model.kappa = bestkappa;
  printf("best parameters: %.2f %.2f %.2f, err/n=%.5f\n", bestalpha, bestbeta, bestkappa, besterror/n);

  blob.model.min = blob.min;
  blob.model.max =blob.max;


  if(error < 1){
    for(float i=blob.position.x; i>0; i-=5){
      if(pdf(blob.position - dvec3(i,0,0)) < 0.001f){
        blob.min.x = i;
        break;
      }
    }
    for(float i=blob.position.x; i<a0; i+=5){
      if(pdf(blob.position + dvec3(i,0,0)) < 0.001f){
        blob.max.x = i;
        break;
      }
    }
    for(float i=blob.position.y; i>0; i-=5){
      if(pdf(blob.position - dvec3(0,i,0)) < 0.001f){
        blob.min.y = i;
        break;
      }
    }
    for(float i=blob.position.y; i<a1; i+=5){
      if(pdf(blob.position + dvec3(0,i,0)) < 0.001f){
        blob.max.y = i;
        break;
      }
    }
    for(float i=blob.position.x; i>0; i-=5){
      if(pdf(blob.position - dvec3(0,0,i)) < 0.001f){
        blob.min.z = i;
        break;
      }
    }
    for(float i=blob.position.x; i<a2; i+=5){
      if(pdf(blob.position + dvec3(0,0,i)) < 0.001f){
        blob.max.z = i;
        break;
      }
    }
  }


  // if (error < 1 ){
  //   while(blob.generalized_multivariate_gaussian_pdf(blob.min) > 0.0001){
  //     blob.min -= vec3(10,10,10);
  //     if(blob.min.x <= 0 && blob.min.y <= 0 && blob.min.z <= 0 )break;
  //   }
  //   while(blob.generalized_multivariate_gaussian_pdf(blob.max) > 0.0001){
  //     blob.max += vec3(10,10,10);
  //     if(blob.max.x >= a0 && blob.max.y >= a1 && blob.max.z >= a2 )break;
  //   }
  // }

  in->model.min = blob.min;
  in->model.max = blob.max;
  in->model.alpha = blob.model.alpha;
  in->model.beta = blob.model.beta;
  in->model.kappa = blob.model.kappa;
  in->model.type = model_type;

  in->model.n = n;
  in->model.error = besterror;

  // blob.min.x = blob.position.x
  return blob;
}

static void kill(ScaleBlob* self){
  self->model.type = '_';
  for(ScaleBlob * c : self->children){
    kill(c);
  }
}

void Estimator::select_scales(Nrrd* source, ScaleBlob* in, ArFilter &filter){

  if(in->n == 0){
    in->model.type = '_';
    for(ScaleBlob *child : in->children){
      // choose the region 
      select_scales(source, child, filter);
    }
    return;
  }
  
  // initialize function
  int a0 = source->axis[0].size;
  int a1 = source->axis[1].size;
  int a2 = source->axis[2].size;

  // choose either self or all children.
  // the children have more detail than the parent.
  // so the region we want to explains is the children.
  // in->model.type = 'g';

  // build a mask containing the pixels
  // that are explained by the children.
  // filter.clear();
  // filter.draw_blobs(in->children, "lm");
  // Nrrd* mask = filter.commit();

  // std::set<ivec3> self_ixs;
  // std::set<ivec3> children_ixs;

  ivec3 min;
  ivec3 max(a0,a1,a2);

  // get the rectangular region which is
  // explained by the children.
  for(ScaleBlob* child : in->children){
    min.x = fmin(min.x, child->min.x);
    min.y = fmin(min.y, child->min.y);
    min.z = fmin(min.z, child->min.z);

    max.x = fmax(max.x, child->max.x);
    max.y = fmax(max.y, child->max.y);
    max.z = fmax(max.z, child->max.z);
  }
  if(min.x<0)min.x=0;
  if(min.y<0)min.y=0;
  if(min.z<0)min.z=0;
  if(max.x>=a0)max.x=a0-1;
  if(max.y>=a1)max.y=a1-1;
  if(max.z>=a2)max.z=a2-1;


  float self_error = 0;
  float children_error = 0;
  // now compute the errors of the parent and children
  // when explaining this region of the data.
  for (int xx=min.x; xx < max.x; xx++){
    for(int yy=min.y; yy < max.y; yy++){
      for(int zz=min.z; zz < max.z; zz++){
        int i = xx + a0*yy + a0*a1*zz;
        vec3 p(xx,yy,zz);

        float region = 0;
        for(ScaleBlob* child : in->children){
          region += child->ellipsepdf(p);
        }

        if (region > 0){
          float val_children = 0;
          float val_self     = 0;
          float overlap      = 0;
          for(ScaleBlob *child : in->children){
            float val = child->generalized_multivariate_gaussian_pdf(p);

            // if(val_children > 0.3f && val > 0.3f){
              // overlap += 1;
            // }
            val_children += val;

          }

          val_self = in->generalized_multivariate_gaussian_pdf(p);
          
          float val_source = ((float*)source->data)[i];

          // if(fabs(val_self-val_source) > 0.3f)self_error += 1;
          // if(fabs(val_children-val_source) > 0.3f)children_error += 1;
          // children_error += overlap;
          self_error += sqrt(fabs(val_self-val_source));
          children_error += sqrt(fabs(val_children - val_source));
          // if(val_children>1)children_error += 1;
        }
      }
    }
  }

  if(self_error > children_error){
    in->model.type = '_';
    for(ScaleBlob *child : in->children){
      // choose the region 
      select_scales(source, child, filter);
    }
  }else{
    kill(in);
    in->model.type = 'g';
  }
  printf(".");
  // for(int i=0;i<)

}