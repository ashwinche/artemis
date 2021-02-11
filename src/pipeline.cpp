#include "pipeline.h"
#include "util.h"
#include <iostream>
#include <glm/glm.hpp>
#include <stdio.h>
#include <string.h>
#include <set>
#include <map>
#include <unordered_map>
#include <queue>
#include <limits>
#include <queue>
#include <SFML/Graphics.hpp>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "blobmath.h"
#include "estimator.h"

#include <unistd.h>

using namespace std::chrono; 
using namespace std;

ReprMode::ReprMode(const char *name){
  this->name = name;
  this->geom = "none";
  scan.v     = 0;
  blob.scale = 0;
  timestep   = 0;
  highlight.blobs = std::vector<ScaleBlob*>();
  highlight.lines = std::vector<vec3>();
  highlight.paths = std::vector<std::vector<ScaleBlob*>>();
  highlight.timestep = -1;
  highlight.locus    = vec3(0,0,0);
  highlight.path_smooth_alpha = 0.9f;
}
bool ReprMode::operator==(ReprMode &r){
  return (name == r.name) && (timestep == r.timestep)
    && ((blob.scale == r.blob.scale));
}
ArPipeline::ArPipeline(ArExperiment *exp){
  this->exp = exp;
  for(int i=0;i<exp->high-exp->low+1;i++){
    ArFrameData data;
    data.blob     = 0;
    data.complete = false;
    // data.tgmm_blob = 0;
    frames.push_back(data);
  }
  store.nbuf = 1;
  store.buf  = new Nrrd*[store.nbuf];
  for(int i=0;i<store.nbuf;++i){
    store.buf[i] = exp->copy(exp->low);
  }
  store.init = true;

  filter.init(store.buf[0]);

  cells = 0;
}
ArFrameData &ArPipeline::get(int frame){
  int i = frame - exp->low;
  while(i >= frames.size() ){
    ArFrameData afd;
    afd.complete = false;
    frames.push_back(afd);
  }
  return frames[i];
}
int ArPipeline::low(){
  return exp->low;
}
int ArPipeline::high(){
  return exp->high;
}

static std::vector<float> collect_scales(ScaleBlob *blob){
  std::set<float> scales;
  std::queue<ScaleBlob*> queue;
  queue.push(blob);
  while(!queue.empty()){
    blob = queue.front();
    queue.pop();
    // printf("blob %p %.2f\n", blob, 14.f);
    // printf("blob %p %.2f\n", blob, blob->scale);
    scales.insert(blob->scale);
    // printf("inserted.\n");
    for(ScaleBlob *child : blob->children){
      queue.push(child);
    }
    // printf("hello\n");
  }
  scales.erase(-1);
  return std::vector<float>(scales.begin(), scales.end());
}
static float compute_epsilon(std::vector<float> &v){
  if(v.size() < 2)return 0;
  float eps = fabs(v[1] - v[0]);
  for(int i=1;i<v.size();++i){
    eps = fmin(eps, fabs(v[i]-v[i-1]));
  }
  return eps/2.f;
}
void ArPipeline::path_highlight(ReprMode *rm, vec3 p, vec3 ray, bool diagnose, bool add){
  if(!get(rm->timestep).complete){
    printf("timestep not yet complete. break.\n");
  }
  ray = normalize(ray)*3.f;
  std::vector<ScaleBlob*> found;
  ArFrameData frame = get(rm->timestep);
  BSPTree<ScaleBlob> *bsptree = &(frame.bspblobs);

  // shoot ray out until it hits a blob.
  // imagine this ray as actually being a "cone" projecting outward,
  // with higher precision near the start.
  for(int i=0;i<100 && found.size()==0;++i){
    p+=ray;
    bsptree->find_within_distance(found, p, i+1);
    if(found.size()>0){
      break;
    }
  }
  ScaleBlob *locus = found[0];
  rm->highlight.highlight_loci.clear();
  rm->highlight.highlight_loci.push_back(found[0]);
}

void ArPipeline::track(ReprMode &repr, std::string output){
  repr.name = "detection";
  // return;
  printf("track.\n");
  repr.geom = "paths";
  repr.name = "cells";

  cells = new std::vector<Cell*>[high() + 1];

  // unlink all blobs:
  for(int i=low();i<=high();i++){
    while(i >= frames.size()){
      ArFrameData data{};
      data.complete = false;
      frames.push_back(data);
    }
    for(ScaleBlob *blob : frames[i].bspblobs.as_vector()){
      blob->succ.clear();
      blob->pred.clear();
    }
  }

  // iterate through frames.
  for(int i=low();i<=high();i++){
    ArFrameData &frame =  get(i);
    if(!frame.complete){
      // printf("load frame.\n");
      loadframe(i);
    }
    if(!frame.complete){
      process(i, i);
      // printf("save frame.\n");
      saveframe(i);
      frame =  get(i);
    }
    std::vector<ScaleBlob*> blobs = frame.bspblobs.as_vector();

    if(i == low()){   // first frame
      for(ScaleBlob* blob : blobs){
        if(blob->scale == frame.scales[0]){
          Cell *cell = new Cell;
          cell->time = 0;
          cell->pred = 0;
          cell->blob = blob;
          cell->f_id = cells[i].size();
          cells[i].push_back(cell);
        }
      }
    }
    else{             // succeeding frames
      std::map<ScaleBlob*, bool> alive;
      for(ScaleBlob *blob : blobs)alive[blob] = true;

      int countdeaths = 0;
      std::unordered_map<Cell*, bool> satisfied; 
      for(Cell* cell : cells[i-1])satisfied[cell] = false;

      printf("frame.\n");
      for(;;){

        // choose the cell to satisfy next, and the blob it connects to.
        Cell *cell = 0;
        ScaleBlob *closest = 0;
        float mindist = 1000.f;

        // for each cell in the previous frame...
        // find the appropriate successor.
        for(Cell *cellc : cells[i-1]){
          if(!cellc)continue;
          if(!cellc->blob)continue;
          if(satisfied[cellc])continue;
          std::vector<ScaleBlob*> potential;
          potential.push_back(frame.bspblobs.find_closest(cellc->blob->position));
          frame.bspblobs.find_within_distance(potential, cellc->blob->position, mindist);
          for(ScaleBlob *blob : potential){
            if(!blob)continue;
            // float dist = glm::distance(cellc->blob->position, blob->position);
            float dist = cellc->blob->wasserstein_distance(blob);
            float rad1 = blob->covmaxev();
            float rad2 = blob->covmaxev();
            // dist -= rad1 + rad2;
            if(dist < mindist){
              cell = cellc;
              mindist = dist;
              closest = blob;
            }
          }
        }
        printf("  d=%.2f\n", mindist);
        // printf("!\n");
        if(!cell)break;
        satisfied[cell] = true;

        


        // if(!closest){
        //   printf("blob died...potential.size = %d\n", potential.size());
        // }
        // if(mindist>2)closest = 0;

        // closest = frame.bspblobs.find_closest(cell->blob->position);

        if(!closest){
          if(cell->life > 0){
            closest = cell->blob;
            Cell *next = new Cell;
            next->f_id = cells[i].size();
            next->pred = cell;
            next->time = i;
            next->blob = closest;
            next->life = cell->life - 1;
            cells[i].push_back(next);
            if(next->pred && next->pred->blob == next->blob){
              printf(".1.ohno!!!");
            }
          }
          else{
            ++countdeaths;
          }
        }
        else{
          Cell *next = new Cell;
          next->f_id = cells[i].size();
          next->pred = cell;
          next->time = i;
          next->blob = closest;
          cells[i].push_back(next);

          cell->blob->succ.push_back(next->blob);
          next->blob->pred.push_back(cell->blob);

          // kill all children/grandchildren
          std::queue<ScaleBlob*> tokill;
          tokill.push(closest);
          while(tokill.size() > 0){
            ScaleBlob *curr = tokill.front();
            tokill.pop();
            alive[curr] = false;
            for(ScaleBlob *child : curr->children){
              tokill.push(child);
            }
          }

          // kill all parents/grandparents
          ScaleBlob *blob = closest;
          while(blob){
            alive[blob] = false;
            blob = blob->parent;
          }
          alive[closest] = false;
          if(next->pred && next->pred->blob == next->blob){
            printf(".2.ohno!!!");
          }
        }
      }
      int countalive = 0;
      for(ScaleBlob *blob : blobs){
        if(alive[blob] && blob->scale == frame.scales[0]){
          printf("cell born ... ");
          ScaleBlob *curr = blob;
          bool valid = true;

          // check that all grandparents are valid.
          while(curr){valid &= alive[curr]; /*printf("%d ", alive[curr]);*/ curr = curr->parent;}
          // printf(" = %d \n", valid);
          if(valid){
            Cell *born = new Cell;
            born->f_id = cells[i].size();

            // find the closest cell.
            for(Cell *cell : cells[i]){

            }
            // ScaleBlob *father = get(i-1).bspblobs.find_closest(blob->position);
            // float dist = glm::distance(cell->blob->position, blob->position);
            // float rad1 = fabs(cbrt(blob->detCov));
            // float rad2 = fabs(cbrt(cell->blob->detCov));
            // dist -= rad1 + rad2;
            // born->pred = dist<0 ? father : 0

            born->newcell = true;
            born->pred = 0;
            born->time = i;
            born->blob = blob;
            cells[i].push_back(born);
            ++countalive;
          }
        }
      }
      printf("%d died. %d born.\n", countdeaths, countalive);
    }


    printf("frame %d has %d cells.\n", i, cells[i].size());
    // write to file


    char *filepath = new char[output.size()+1];
    strcpy(filepath, output.c_str());
    // char filepath[] = "/home/ashwin/data2/tgmmout/17-05-01-4-?????.xml";    // output file format
    char *iwildcard  = filepath;                                          // pointer to wildcard ?s
    int sizewildcard = 5;                                                 // number of ?s
    while(*iwildcard != '?')++iwildcard;

    sprintf(iwildcard, "%05d.xml", i);
    printf("write to file %s\n", filepath);
    FILE *file = fopen(filepath, "w");

    if(file){

      fprintf(file, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<document>");

      for(Cell *cell : cells[i]){
        ScaleBlob *blob = cell->blob;
        int m_id = cell->f_id;
        int m_parent = cell->pred? cell->pred->f_id : -1;

        float m_x = blob->position.x;
        float m_y = blob->position.y;
        float m_z = blob->position.z;

        fprintf(file, 
          "<GaussianMixtureModel id=\"%d\" parent=\"%d\" m=\"%f %f %f\">\n</GaussianMixtureModel>\n",
           m_id, m_parent, m_x, m_y, m_z);
      }

      fprintf(file, "</document>");
      fclose(file);
    }
    delete[] filepath;

  }

  // exit(0);

  // for(int i=low();i<high();i++){
  //   printf("capture %d\n", i);
  //   filter.capture(exp->get(i));

  //   DiscreteKernel kernel = filter.gaussian(2, 2*4);
  //   filter.set_kernel(kernel);
  //   filter.filter();
  //   kernel.destroy();
  //   // filter.highlight(filter.find_maxima());
  //   // filter.commit(store.buf[0]);
  //   // return;
  //   frames.push_back(Frame());
  //   printf("find blobs.\n");
  //   frames[i].blobs    = filter.find_blobs();
  //   frames[i].bspblobs = filter.get_bsp(10);
  //   for(ScaleBlob *blob : frames[i].blobs){
  //     frames[i].bspblobs.insert(blob, blob->mode);
  //   }
  //   printf("connect.\n");
  //   // connect to the previous frame.
  //   if(i>1){
  //     for(ScaleBlob *succ : frames[i].blobs){
  //       ScaleBlob *pred = frames[i-1].bspblobs.find_closest(succ->mode);
  //       if(pred->distance(succ) < 2.f){
  //         pred->succ.push_back(succ);
  //         succ->pred.push_back(pred);
  //       }
  //     }
  //   }

  //   printf("output.\n");
  //   Frame &f = frames[i];



  //   fprintf(file, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<document>");
    

  //   // get IDs for current frame.
  //   pointer_to_id[currmap].clear();
  //   for(int j=0;j<f.blobs.size();j++){
  //     pointer_to_id[currmap][f.blobs[j]] = j;
  //   }


  //   for(ScaleBlob *blob : f.blobs){
  //     int m_id = pointer_to_id[currmap][blob];
  //     int m_parent = -1;
  //     if(blob->pred.size() > 0){
  //       m_parent = pointer_to_id[!currmap][blob->pred[0]];
  //     }
  //     float m_x = blob->position.x;
  //     float m_y = blob->position.y;
  //     float m_z = blob->position.z;
  //     fprintf(file, 
  //       "<GaussianMixtureModel id=\"%d\" parent = \"%d\" m = \"%f %f %f\">\n</GaussianMixtureModel>\n",
  //        m_id, m_parent, m_x, m_y, m_z);
  //   }

  //   fprintf(file, "</document>");
  //   fclose(file);

  //   currmap = !currmap;
  // }

  // for(Frame &f : frames){
  //   for(ScaleBlob *blob : f.blobs){
  //     delete blob;
  //   }
  // }
}

void ArPipeline::repr_highlight(ReprMode *rm, vec3 p, vec3 ray, bool drawpaths, bool add){

  if(get(rm->timestep).complete){
    // printf("highlight along (%.2f %.2f %.2f) + (%.2f %.2f %.2f)\n",p.x,p.y,p.z, ray.x,ray.y,ray.z);
    ray = normalize(ray)*3.f;
    std::vector<ScaleBlob*> found;
    ArFrameData frame = get(rm->timestep);
    BSPTree<ScaleBlob> *bsptree = &(frame.bspblobs);

    printf("number of blobs = %d\n", bsptree->n);

    // shoot ray out until it hits a blob.
    for(int i=0;i<100 && found.size()==0;++i){
      p+=ray; 
      bsptree->find_within_distance(found, p, i+1);
      if(found.size()>0){
        // if(drawpaths){
        //   printf("highlight blob %p:\n", found[0]);
        //   found[0]->print();
        // }
        if(!add){
          rm->highlight.highlight_loci.clear();
        }
        rm->highlight.highlight_loci.push_back(found[0]);

        found = std::vector<ScaleBlob*>();
        // found = bsptree->as_vector();
        bsptree->find_within_distance(found, p, 625);
      }
    }

    rm->highlight.blobs.clear();
    // for(int i=0;i<found.size();i++)rm->highlight.blobs.push_back(found[i]);

    bool drawblobs = true;

    if(drawpaths){
      printf("drawpaths %.2f %.2f %.2f.\n", p.x, p.y, p.z);
      if(add){
        // rm->highlight.paths = paths;
      }else{
        rm->highlight.paths.clear();
      }
      bool found = false;
      for(std::vector<ScaleBlob*> path : paths){
        // if(found)break;
        // if(path[0] == found[0]){
        //   rm->highlight.paths.push_back(path);
        // }
        for(ScaleBlob *blob : path){
          if(!blob)break;
          // printf("blob=%p\n,", blob);
          if(glm::distance(vec3(blob->position), p) < 10){
            rm->highlight.paths.push_back(path);
            if(drawblobs){
              for(ScaleBlob *blob2 : path){
                rm->highlight.blobs.push_back(blob2);
              }
              // found = true;
              // drawblobs = false;
            }
            break;
          }
          // break;
        }
      }

    }else{
      rm->highlight.paths.clear();
    }
    printf("draw blobs %d\n", rm->highlight.blobs.size());
    rm->highlight.locus = p;
    rm->highlight.timestep = rm->timestep;
    // rm->highlight.paths = longest_paths(rm->highlight.blobs);
    // printf("highlight %d\n",found.size());
  }
  // rm->highlight.lines.push_back(p+ray*1.f);
  // rm->highlight.lines.push_back(p+ray*100.f);
}

ReprMode ArPipeline::repr_coarser(ReprMode rm){
  std::vector<float> scales = get(rm.timestep).scales;
  ++rm.blob.scale;
  if(rm.blob.scale < 0)rm.blob.scale = 0;
  if(rm.blob.scale >= scales.size()) rm.blob.scale = scales.size()-1;
  return rm;
}
ReprMode ArPipeline::repr_finer(ReprMode rm){
  std::vector<float> scales = get(rm.timestep).scales;
  --rm.blob.scale;
  if(rm.blob.scale < 0)rm.blob.scale = 0;
  if(rm.blob.scale >= scales.size()) rm.blob.scale = scales.size()-1;
  return rm;
}

static std::vector<ScaleBlob*> collect_blobs(ScaleBlob *blob, float scalemin, float scalemax){
  // scalemin = -1000;
  // scalemax = 10000;collect_blobs
  // printf("(");
  // printf("collect_blobs %.2f %.2f\n", scalemin, scalemax);
  std::vector<ScaleBlob*> blobs;
  // printf("blob = %p\n", blob);
  // printf("blob.scale = %.2f\n", blob->scale);
  if(blob->scale >= scalemin || blob->scale <= 0 || scalemin == -1){
    if(blob->scale <= scalemax || blob->scale <= 0 || scalemax == -1){
      // printf(".");
      // printf("blob.scale = %.2f\n", blob->scale);
      blobs.push_back(blob);
    }
    // printf("  bscale = %.2f, %.2f\n", blob->scale, scalemax);
    for(ScaleBlob *child : blob->children){
      std::vector<ScaleBlob*> childblobs = collect_blobs(child, scalemin, scalemax);
      blobs.insert(blobs.end(), childblobs.begin(), childblobs.end());
    }
  }
  // printf(")");
  return blobs;
}

static std::vector<ScaleBlob*> modeled_blobs(ScaleBlob *blob){
  std::vector<ScaleBlob*> blobs;
  if(blob->model.type != '_'){
    blobs.push_back(blob);
  }
  for(ScaleBlob *child : blob->children){
    std::vector<ScaleBlob*> childblobs = modeled_blobs(child);
    blobs.insert(blobs.end(), childblobs.begin(), childblobs.end());
  }

  return blobs;
}


static std::vector<ScaleBlob*> detected_blobs_ch(ScaleBlob *self){
  bool add_self = true;
  std::vector<ScaleBlob*> childblobs;
  float children_n;
  float children_error = 2.f;

  // add up blobs from all of the children.
  for(ScaleBlob* child : self->children){
    std::vector<ScaleBlob*> next = detected_blobs_ch(child);
    childblobs.insert(childblobs.end(), next.begin(), next.end());
    for(ScaleBlob* blob : next){
      children_error = fmin(children_error, blob->model.error/blob->model.n);
    }
  }

  // check if the children have a better model than we do.
  if(children_error < self->model.error / self->model.n){
    add_self = false;
  }

  if(self->n == 0)add_self = false;

  // if not, then just return our model. 
  if(add_self || self->children.size() == 0){
    std::vector<ScaleBlob*> selfblob;
    selfblob.push_back(self);
    // printf("self ");
    return selfblob;
  }
  else{
    // printf("chrn ");
    // otherwise, return the childrens' model.
    return childblobs;
  }
}
static std::vector<ScaleBlob*> detected_blobs(ScaleBlob *blob){
  std::vector<ScaleBlob*> ret = detected_blobs_ch(blob);
  // printf("%d\n", ret.size());
  return ret;
  std::vector<ScaleBlob*> blobs;

  bool add_children = false;
  if(blob->n == 0)add_children = true;

  if(blob->children.size() > 0){
    float error_this = blob->model.error/blob->model.n;
    float sum_err_children = 0;
    float sum_n_children = 0;
    for(ScaleBlob *child : blob->children){
      sum_err_children += child->model.error;
      sum_n_children   += child->model.n;
    }
    float error_children = sum_err_children/sum_n_children;

    // printf("errors = %.4f %.4f\n", error_children, error_this);

    if(error_children < error_this){
      add_children = true;
    }
  }
  if(add_children){
    for(ScaleBlob *child : blob->children){
      std::vector<ScaleBlob*> childblobs = detected_blobs(child);
      blobs.insert(blobs.end(), childblobs.begin(), childblobs.end());
    }
  }else{
    blobs.push_back(blob);
  }
  // printf(")");
  return blobs;
}

/* find all paths in the part of the dataset that has been
 * processed so far, so long as the path is greater than
 * minlen in length.
 *
 */

static bool fun_sort_blob_by_n(ScaleBlob* a, ScaleBlob* b){
  return a->n > b->n;
}
void ArPipeline::findpaths(int minlen, int maxframe, const char* mode){
  if(maxframe <= 0)maxframe = frames.size();
  if(maxframe > frames.size())maxframe = frames.size();
  std::vector<ScaleBlob*> allblobs;
  for(int i=0;i<maxframe;i++){   // for each frame
    if(!strcmp(mode, "tgmm") || frames[i].complete){           // if it is processed
                                      // collect blobs
      std::vector<ScaleBlob*> blobsi;
      blobsi = collect_blobs(frames[i].blob, 0, std::numeric_limits<float>::infinity());
      std::sort(blobsi.begin(), blobsi.end(), fun_sort_blob_by_n);
      // for(auto sb : blobsi){
      //   printf("%.2f; ", sb->n);
      // }printf("\n");
      allblobs.insert(allblobs.end(), blobsi.begin(), blobsi.end());
      // minlen = int(i*0.75);
    }
  }
  printf("find paths > %d.\n", minlen);
  printf("find paths in %lu blobs.\n", allblobs.size());
  std::vector<std::vector<ScaleBlob*>> allpaths;
  if(!strcmp(mode, "quick"))   allpaths = longest_paths2(allblobs, minlen);
  if(!strcmp(mode, "tgmm"))    allpaths = longest_paths2(allblobs, minlen);
  if(!strcmp(mode, "longest")) allpaths = longest_paths2(allblobs, minlen);
  paths.clear();
  for(int i=0;i<allpaths.size();++i){
    if(allpaths[i].size() > minlen){
      paths.push_back(allpaths[i]);
    }
  }
  paths = allpaths;
  printf("found %lu paths.\n", paths.size());
}

// attempt to parallelize it. doesn't work though.
    // static void process_frame(ArExperiment *exp, ArFrameData *frame, int framen){
    //   ArFilter filter;
    //   filter.init(exp->get(framen));
    //   // filter.capture(exp->get(framen));
    //   printf("process %d\n", framen);
    //   ScaleBlob *blob             = filter.compute_blob_tree();


    //   std::vector<float> scales   = collect_scales(blob);
    //   frame->blob      = blob;
    //   frame->scales    = scales;
    //   frame->scale_eps = compute_epsilon(scales);
    //   frame->complete  = true;
    //   frame->bspblobs  = filter.get_bsp(10);

    //   BSPTree<ScaleBlob> *bsptree = &(frame->bspblobs);
    //   std::vector<ScaleBlob*> allblobs = collect_blobs(blob, 0, std::numeric_limits<float>::infinity());
      
    //   for(ScaleBlob *sb : allblobs){
    //     bsptree->insert(sb, sb->position);
    //   }
    //   filter.destroy();
    // }

// void ArPipeline::link_frames(int low, int high){
//   if(high > exp->high){
//     high = exp->high;
//     printf("pipeline.link_frames, low=%d, high=%d\n", low, high);
//   }
//   for(int frame=low;frame<=high-1;frame++){
//     if(!get(frame).complete)break;
    
//     BSPTree<ScaleBlob> *t0 = &frames[ frame - exp->low    ].bspblobs;
//     BSPTree<ScaleBlob> *t1 = &frames[ frame - exp->low +1 ].bspblobs;

//     std::vector<ScaleBlob*> v0 = t0->as_vector();
//     for(ScaleBlob *sb : v0){
//       std::vector<ScaleBlob*> potential;
//       t1->find_within_distance(potential, sb->position, 1000.f);
//       for(int i=0;i<potential.size();i++){

//         //  *** NOTE ***
//         // 1.0 means touching. I don't know what 2.0 means. This is somewhat
//         // arbitrary, but we notice that in the data, sometimes successor
//         // blobs aren't touching predecessor blobs, so we need a small margin.
//         // In any case, we next want to choose the potential successor with
//         // the smallest distance and closest scale.
//         //  *** **** ***

//         if(sb->n>1 && potential[i]->n>1 && sb->distance(potential[i]) <= 1.f){
//           if(sb->detCov/potential[i]->detCov < 2.f && potential[i]->detCov/sb->detCov < 8.f){
//             // volume cannot more than double or half.
//             sb->succ.push_back(potential[i]);
//             potential[i]->pred.push_back(sb);
//           }
//         }
//         // printf("%d.", i);
//       }
//       // printf("o");
//       ++itr;
//       // for(ScaleBlob *sb0 : sb->succ){
//       //   sb0->pred.push_back(sb);
//       // }
//     }
//   }
// }
/* store.buf[0] contains the output of the filter chain.
 * store.buf[0] contains the output of repr().
 */

void ArPipeline::link(int low, int high){
  printf("linking...");
  for(int frame = this->low(); frame <= high; ++frame){
    printf("%d ", frame);
    if(get(frame).complete && get(frame+1).complete);
    else continue;
    BSPTree<ScaleBlob> *t0 = &frames[ frame - exp->low    ].bspblobs;
    BSPTree<ScaleBlob> *t1 = &frames[ frame - exp->low +1 ].bspblobs;

    std::vector<ScaleBlob*> v0 = t0->as_vector();
    std::vector<ScaleBlob*> v1 = t1->as_vector();
    // clear predecessors so we can rewrite them.
    for(ScaleBlob *sb : v1){
      sb->pred.clear();
    }
    for(ScaleBlob *sb : v0){
      sb->succ.clear();       // clear successors so we can rewrite them.
      std::vector<ScaleBlob*> potential;
      ScaleBlob* closest = t1->find_closest(sb->position);
      potential.push_back(closest);
      t1->find_within_distance(potential, sb->position, 1000.f);
      for(int i=0;i<potential.size();i++){
        if(sb->n>1 && potential[i]->n>1 && sb->distance(potential[i]) <= 1.f){
          // if(sb->detCov/potential[i]->detCov < 8.f && potential[i]->detCov/sb->detCov < 8.f){            // volume cannot more than double or half.
          sb->succ.push_back(potential[i]);
          potential[i]->pred.push_back(sb);
          // }
        }
      }
      std::sort(sb->succ.begin(), sb->succ.end(), fun_sort_blob_by_n);
    }
  }
  printf("done.\n");

}
void ArPipeline::process(int low, int high){
  if(high > exp->high){
    high = exp->high;
    printf("pipeline.process, low=%d, high=%d\n", low, high);
  }
  for(int frame=low;frame<=high;frame++){
    printf("pipeline.process %d\n",frame);
    tick("");

    // THIS CAN BE PARALLELIZED {

        // ARFilter filter;
        filter.capture(exp->get(frame));
        // printf("process %d\n", frame)
        ScaleBlob *blob             = filter.compute_blob_tree();


        std::vector<float> scales   = collect_scales(blob);
        frames[frame - exp->low].blob      = blob;
        frames[frame - exp->low].scales    = scales;
        frames[frame - exp->low].scale_eps = compute_epsilon(scales);
        frames[frame - exp->low].complete  = true;
        frames[frame - exp->low].bspblobs  = filter.get_bsp(10);

        BSPTree<ScaleBlob> *bsptree = &frames[frame - exp->low].bspblobs;
        std::vector<ScaleBlob*> allblobs = collect_blobs(blob, -1, std::numeric_limits<float>::infinity());
        
        for(ScaleBlob *sb : allblobs){
          bsptree->insert(sb, sb->mode);
        }

    // } THIS CAN BE PARALLELIZED

    // link with previous frame.
    // printf("linking.\n");
    // if(frame-1 >= this->low() && get(frame-1).complete){
    //   BSPTree<ScaleBlob> *t0 = &frames[ frame - exp->low -1 ].bspblobs;
    //   BSPTree<ScaleBlob> *t1 = &frames[ frame - exp->low    ].bspblobs;

    //   std::vector<ScaleBlob*> v0 = t0->as_vector();
    //   int itr = 0;
    //   for(ScaleBlob *sb : v0){
    //     std::vector<ScaleBlob*> potential;
    //     t1->find_within_distance(potential, sb->position, 1000.f);
    //     for(int i=0;i<potential.size();i++){

    //       //  *** NOTE ***
    //       // 1.0 means touching. I don't know what 2.0 means. This is somewhat
    //       // arbitrary, but we notice that in the data, sometimes successor
    //       // blobs aren't touching predecessor blobs, so we need a small margin.
    //       // In any case, we next want to choose the potential successor with
    //       // the smallest distance and closest scale.
    //       //  *** **** ***

    //       if(sb->n>1 && potential[i]->n>1 && sb->distance(potential[i]) <= 1.f){
    //         if(sb->n/potential[i]->n < 8.f && potential[i]->n/sb->n < 8.f){
    //           // volume cannot more than double or half.
    //           sb->succ.push_back(potential[i]);
    //           potential[i]->pred.push_back(sb);
    //         }
    //       }
    //       // printf("%d.", i);
    //     }
    //     std::sort(sb->succ.begin(), sb->succ.end(), fun_sort_blob_by_n);

    //     // for(ScaleBlob *next : sb->succ){
    //     //   printf("%.2f; ", next->n);
    //     // }printf("\n");
    //     // printf("o");
    //     ++itr;
    //     // for(ScaleBlob *sb0 : sb->succ){
    //     //   sb0->pred.push_back(sb);
    //     // }
    //   }
    // }
    tick("done.\n");
    saveframe(frame);
  }
  filter.commit(store.buf[0]);
}

ArGeometry3D* ArPipeline::reprgeometry(ReprMode &mode){
  geometry.lines = std::vector<vec3>();
  geometry.lines_c = std::vector<sf::Color>();
  // printf("mode %s\n", mode.geom);
  if(!strcmp(mode.geom, "paths")){
    // draw reference x- y- z- axes 
    geometry.lines.push_back(dvec3(0,0,0));
    geometry.lines_c.push_back(sf::Color(255,0,0,255));
    geometry.lines.push_back(dvec3(100,0,0));
    geometry.lines_c.push_back(sf::Color(255,0,0,255));
    geometry.lines.push_back(dvec3(0,0,0));
    geometry.lines_c.push_back(sf::Color(0,255,0,255));
    geometry.lines.push_back(dvec3(0,100,0));
    geometry.lines_c.push_back(sf::Color(0,255,0,255));
    geometry.lines.push_back(dvec3(0,0,0));
    geometry.lines_c.push_back(sf::Color(0,0,255,255));
    geometry.lines.push_back(dvec3(0,0,100));
    geometry.lines_c.push_back(sf::Color(0,0,255,255));
    
    // draw trajectory of highlighted blobs
    if(!mode.highlight.paths.empty()){

      // draw the longest trajectories.
      // double path_smooth_beta = 1.f - mode.highlight.path_smooth_alpha;
      // printf("Draw %d paths.\n", mode.highlight.paths.size());

      for(std::vector<ScaleBlob*> path : mode.highlight.paths){

        std::vector<glm::dvec3> smoothed(path.size());

        if(true){   // smooth path
          int smooth_path = 0;
          for(int i=0;i<path.size();i++){
            int k = smooth_path;
            if(i-k < 0){
              k = i;
            }
            if(i+k >= path.size()){
              k = path.size() - i - 1;
            }
            int j0 = i-k;
            int j1 = i+k;
            float n=k*2+1;
            for(int j=j0;j<=j1;++j){
              smoothed[i] += path[j]->position;
            }
            smoothed[i] /= n;
            // smoothed[i] = path[i]->position;
          }
        }else{
          // smoothed = path;
        }

        if(false){  // draw vectors
          smoothed = std::vector<glm::dvec3>(2);
          smoothed[0] = path[0]->position;
          smoothed[1] = path[path.size()-1]->position;
          // smoothed[1] = path[path.size()-1]->position - path[1]->position;
          // smoothed[1] *= 20.f/glm::length(smoothed[1]);
          // smoothed[1] += smoothed[0];
        }

        float len  = float(smoothed.size());
        // if(len<20)continue;
        float step = 1.f/len;
        glm::dvec3 weightedp = path[0]->position;
        for(int j=0;j<smoothed.size()-1;j++){
          // geometry.lines.push_back(path[j]->position);
          // geometry.lines.push_back(path[j+1]->position);
          // geometry.lines.push_back(weightedp);
          // glm::dvec3 weightedq = path[j+1]->position;
          
          // weightedp = (weightedp * mode.highlight.path_smooth_alpha) + (weightedq * path_smooth_beta);

          // geometry.lines.push_back(weightedp);

          geometry.lines.push_back(smoothed[j]);
          geometry.lines.push_back(smoothed[j+1]);

          const int SMOOTH   = 0;
          const int GRADIENT = 1;
          int pathcolormode = GRADIENT;

          if(pathcolormode == SMOOTH){
            int r0 = 0   + int(200.f * step * j    );
            int g0 = 155 + int(100.f * step * j    );
            int b0 = 255 - int(250.f * step * j    );
            int r1 = 0   + int(200.f * step * (j+1));
            int g1 = 155 + int(100.f * step * (j+1));
            int b1 = 255 - int(250.f * step * j    );
            geometry.lines_c.push_back(sf::Color(r0,g0,b0,255));
            geometry.lines_c.push_back(sf::Color(r1,g1,b1,255));
          }
          else if(pathcolormode == GRADIENT){
            glm::dvec3 direction = path[path.size()-1]->position - path[0]->position;
            float len = glm::length(direction);
            direction = glm::normalize(direction);
            int alpha = len*15;
            if(alpha>255)alpha = 255;
            int r,g,b;
            r = 127 + direction.x*127;
            g = 127 + direction.y*127;
            b = 127 + direction.z*127;
            geometry.lines_c.push_back(sf::Color(r,g,b,alpha));
            geometry.lines_c.push_back(sf::Color(r,g,b,alpha));
          }
        }
      }
    }
  }
  if(cells){
    // geometry.lines.push_back(dvec3(0,0,0));
    // geometry.lines_c.push_back(sf::Color(255,0,0,255));
    // geometry.lines.push_back(dvec3(100,100,100));
    // geometry.lines_c.push_back(sf::Color(255,0,0,255));
    // geometry.lines.push_back(dvec3(0,0,0));
    // geometry.lines_c.push_back(sf::Color(0,255,0,255));
    // geometry.lines.push_back(dvec3(0,100,0));
    // geometry.lines_c.push_back(sf::Color(0,255,0,255));
    // geometry.lines.push_back(dvec3(0,0,0));
    // geometry.lines_c.push_back(sf::Color(0,0,255,255));
    // geometry.lines.push_back(dvec3(0,0,100));
    // geometry.lines_c.push_back(sf::Color(0,0,255,255));
    // draw predecessors of cells

    ArFrameData frame = get(mode.timestep);
    std::vector<Cell*> cellsi = cells[int(mode.timestep)];
    // printf("cellsi: %d\n", cellsi.size());

    std::vector<ScaleBlob*> draw;
    std::vector<ScaleBlob*> color;
    for(Cell *cell : cellsi){
      // printf(".");
      // singleton[0] = cell->blob;
      if(cell->pred){
        // printf("pred!\n");
        geometry.lines.push_back(cell->pred->blob->position);
        geometry.lines.push_back(cell->blob->position);
        geometry.lines_c.push_back(sf::Color(255,50,255,255));
        geometry.lines_c.push_back(sf::Color(255,50,255,255));
        // if(glm::distance(cell->pred->blob->position,cell->blob->position) < 0.000001){
        //   // printf("%p %p\n", cell, cell->pred);
        //   // printf("d = %.5f\n", glm::distance(cell->pred->blob->position,cell->blob->position));
        //   // printf("d = %.2f %.2f %.2f ; %.2f %.2f %.2f\n", cell->pred->blob->position.x,cell->pred->blob->position.y,cell->pred->blob->position.z,cell->blob->position.x,cell->blob->position.y,cell->blob->position.z);
        //   geometry.lines.push_back(dvec3(0,0,0));
        //   geometry.lines.push_back(cell->blob->position);
        //   geometry.lines_c.push_back(sf::Color(255,200,100,255));
        //   geometry.lines_c.push_back(sf::Color(255,200,100,255));
        // }
      }
      // printf("\n");
    }
  }

  for(vec3 v : mode.highlight.lines){
    geometry.lines.push_back(v);
    geometry.lines_c.push_back(sf::Color::Yellow);
  }

  return &geometry;
}

void ArPipeline::estimate(int timestep){
  ArFrameData frame = get(timestep);
  if(!frame.complete)return;
  std::vector<ScaleBlob*> blobs = collect_blobs(frame.blob, -1, -1);

  Estimator estimator;

  // filter.capture(exp->get(timestep));

  filter.capture(exp->get(timestep));
  float scale = frame.scales[0];
  filter.lapofgaussian_masked(scale, false);
  filter.commit(store.buf[0]);

  ScaleBlob estimated_blob;
  printf("estimating %lu blobs.\n", blobs.size());
  for(int i=0;i<blobs.size();i++){
    printf("%d. ", i);
    estimated_blob = estimator.fit(store.buf[0], blobs[i]);
    // blobs[i]->model.type = 'g';
  }

  // for(int i=100;i<blobs.size();i++){
  //   blobs[i]->model.kappa = 0;
  // }
  printf("done.\n");
}
void ArPipeline::select_scales(int timestep){
  ArFrameData frame = get(timestep);
  if(!frame.complete)return;
  Estimator estimator;

  filter.capture(exp->get(timestep));
  float scale = frame.scales[0];
  filter.lapofgaussian_masked(scale, false);
  filter.commit(store.buf[0]);

  estimator.select_scales(store.buf[0], frame.blob, filter);
}
/**
 * Chooses how to represent the volume-rendered data.
 * Eg. we may represent the data as:
 *   - "plain"     : raw data
 *   - "blobs"     : depicting the blobs recovered from the data
 *   - "gaussian"  : rendering a blurred version of the data at some sigma
 *   - "laplacian" : laplacian of gaussian filter at some sigma
*/
Nrrd *ArPipeline::repr(ReprMode &mode, bool force){
  if(mode.timestep < exp->low)mode.timestep = exp->low;
  if(mode.timestep > exp->high)mode.timestep = exp->high;
  if(!strcmp(mode.name, "plain")){
    // plain representation is simple and is always allowed.
    // default to plain if there the frame has not been processed yet.
    // printf("repr plain\n");
    return exp->get(mode.timestep);
  }

  if(!strcmp(mode.name, "filter residue")){
    printf("repr %s\n", mode.name);
    return store.buf[2];
  }
  if(!strcmp(mode.name, "filter internal")){
    printf("repr %s\n", mode.name);
    return store.buf[0];
  }
  if(!strcmp(mode.name, "sandbox")){
    int last_np  = 0;
    int last_dnp = 0;
    int last_ddnp = 0;
    float f = 0.5f;
    if(f>6)f = 0.5f;
    else f *= 1.1f;
    // else f += 0.3f;
    // for(float f = 0.5f; f< 6.f; f*= 1.33f){
    tick("sandbox");
    for(float f = 3.6; f< 3.7; f++){
      filter.capture(exp->get(mode.timestep));

      DiscreteKernel kernel = filter.gaussian(f, int(f*4));
      
      filter.max1();    
      filter.set_kernel(kernel);
      filter.filter();
      tick("gaussian");
      kernel.destroy();
      filter.laplacian3d();
      tick("laplacian");
      // filter.find_blobs();
      // filter.threshold(0.01,1.f);
      // filter.filter();
      // int np  = filter.find_maxima().size();
      // std::vector<ScaleBlob*> blobs = filter.find_blobs();
      int np = filter.count_blobs();
      tick("count");
      // int np  = filter.count_connected_components();
      // int np  = 0;
      int dnp = last_np - np; 
      int ddnp = last_dnp - dnp;
      int dddnp = last_ddnp - ddnp;
      // printf("%d\n", np);
      // printf("%.2f\t%d\t%d\t%d\t%d\n", f, np, dnp, ddnp, dddnp);
      printf("%.2f\t%d\n", f, np);
      if(dnp > 0){
        printf(" **** \n");
        break;
      }
      last_np = np;
      last_dnp = dnp;
      last_ddnp = ddnp;
    }
    filter.normalize();
    filter.commit(store.buf[0]);
    return store.buf[0];
  }
  // if(!get(mode.timestep).complete){
  //   return exp->get(mode.timestep);
  // }

  static ReprMode last_repr("");
  static Nrrd *last_nrrd;
  if(mode == last_repr && !force){
    // printf("repr %s unchanged.\n",mode.name);
    return last_nrrd;
  }
  last_repr  = mode;
  // printf("repr %s\n",mode.name);

  if(!strcmp(mode.name, "dconv")){
    // deconvolved image, using richardson-lucy.
    filter.capture(exp->get(mode.timestep));
    filter.deconvolve();
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }

  if(!strcmp(mode.name, "detection")){
    filter.capture(exp->get(mode.timestep));
    DiscreteKernel kernel = filter.gaussian(2, 2*4);
    filter.set_kernel(kernel);
    filter.filter();
    kernel.destroy();
    std::vector<ScaleBlob*> blobs = filter.find_blobs();
    filter.draw_blobs(blobs, ".+");
    filter.commit(store.buf[0]);
    for(ScaleBlob *b : blobs)delete b;
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "sum")){
    filter.clear();
    for(int i=0;i<30;i++){
      printf("%d ", i);
      ArFrameData frame = get(i);
      float scalemin = frame.scales[mode.blob.scale] - frame.scale_eps;
      float scalemax = frame.scales[mode.blob.scale] + frame.scale_eps;
      std::vector<ScaleBlob*> blobs = collect_blobs(frame.blob, scalemin, scalemax);      filter.clear();
      filter.draw_blobs(blobs, "g+");
    }
    printf("\n");
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "flow")){
    // visualize the flow of GFP throughout the experiment.
    
    printf("repr %s\n", mode.name); 
    ArFrameData frame = get(mode.timestep);   
    filter.capture(exp->get(exp->low));

    // use a gaussian kernel
    float scale = frame.scales[mode.blob.scale];
    DiscreteKernel kernel = filter.gaussian(scale, int(scale*4));
    filter.set_kernel(kernel);

    // perform laplacian of gaussian
    tick("flow");
    filter.max1();
    filter.filter();
    filter.laplacian3d();
    filter.normalize();

    // get maxima at t=0
    std::vector<std::vector<glm::vec3>> paths;

    std::vector<ivec3> maxima = filter.find_maxima();
    for(ivec3 maximum : maxima){
      std::vector<vec3> path;
      path.push_back(vec3(maximum.x, maximum.y, maximum.z));
      paths.push_back(path);
    }

    // highlight.push_back(filter.find_maxima());

    tick("t=0");

    for(int i=0;i<10;i++){
      printf("lap of gaus\n");
      // laplacian of gaussian of frame i
      filter.capture(exp->get(exp->low + i));
      filter.max1();          // max
      filter.filter();        // gaussian
      filter.laplacian3d();   // laplacian
      filter.normalize();     // normalize
      printf("hill climb\n");
      // hill-climb into the next frame.
      for(int p=0;p<paths.size();p++){
        ivec3 point(paths[p][i].x, paths[p][i].y, paths[p][i].z);
        ivec3 hill = filter.hill_climb(point);
        paths[p].push_back(hill);
      }
      // tick("hill-climb");
    }

    printf("found %lu paths.\n", paths.size());
    filter.clear();
    for(int i=0;i<paths.size();i++){
      for(int j=0;j<paths[i].size()-1;j++){
        // printf("  line from %.2f %.2f %.2f to %.2f %.2f %.2f\n",
         // paths[i][j].x,   paths[i][j].y,   paths[i][j].z,
         // paths[i][j+1].x, paths[i][j+1].y, paths[i][j+1].z);
         filter.rasterlineadd(paths[i][j], paths[i][j+1], 1.f, 1.f);
      }
    }
    printf("gaussian.\n");

    DiscreteKernel kernel2 = filter.gaussian(1.f, 5);
    filter.set_kernel(kernel2);
    filter.filter();        // gaussian
    filter.normalize();
    kernel2.destroy();

    kernel.destroy();

    // filter.clear();
    // filter.highlight(highlight[0]);
    // filter.highlight(highlight[1]);

    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "blobs_all")){
    

    ArFrameData frame = get(mode.timestep);
    
    // float scale = frame.scales[mode.blob.scale];
    // filter.capture(exp->get(mode.timestep));
    // DiscreteKernel kernel = filter.gaussian(scale/1, int(scale*4));
    // filter.set_kernel(kernel);
    // filter.max1();
    // filter.filter();
    // filter.laplacian3d();
    // std::vector<ScaleBlob*> blobs = filter.find_blobs();
    // filter.clear();
    // kernel.destroy();

    float scalemin = frame.scales[mode.blob.scale] - frame.scale_eps;
    float scalemax = frame.scales[mode.blob.scale] + frame.scale_eps;
    // printf("scale = %d %.2f %.2f\n", mode.blob.scale, scalemin, scalemax);
    std::vector<ScaleBlob*> blobs = collect_blobs(frame.blob, scalemin, scalemax);

    printf("collect %d blobs\n", blobs.size());
    for(ScaleBlob *sb : blobs){
      // printf("%.2f %.2f %.2f %.2f\n", sb->position.x, sb->position.y, sb->position.z, sb->n);
    }

    // filter.clear();

    // return (last_nrrd = store.buf[0]);
    // float scalemin = frame.scales[mode.blob.scale] - frame.scale_eps;
    // float scalemax = frame.scales[mode.blob.scale] + frame.scale_eps;
    // printf("mode %s. view scale %d with %.2f %.2f\n", mode.name, scalemin, scalemax);
    


    // std::vector<ScaleBlob*> blobs = modeled_blobs(frame.blob);
    

    // std::vector<ScaleBlob*> blobs = collect(frame.blob);

    filter.clear();
    // filter.draw_blobs(blobs, "lm");
    filter.draw_blobs(blobs, "M+");
    filter.normalize();
    // filter.difference_image(exp->get(mode.timestep));
    // filter.draw_blobs(blobs, ".m");

    // std::vector<ScaleBlob*> only;
    // Estimator estimator;
    // // filter.clear();
    // ScaleBlob estimated_blob;
    // only.push_back(&estimated_blob);
    // printf("estimating %lu blobs.\n", blobs.size());
    // for(int i=0;i<blobs.size();i++){
    //   estimated_blob = estimator.fit(exp->get(mode.timestep), blobs[i]);
    //   filter.draw_blobs(only, "m+");
    //   filter.draw_blobs(only, ".m");
    // }

    
    
    printf("done blobs_all.\n");
    // filter.difference_image(exp->get(mode.timestep));
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);   
  }
  if(!strcmp(mode.name, "tgmm")){
    // printf("mode = tgmm\n");
    // ArFrameData frame = get(mode.timestep);
    // // printf("size[%d] = %p = 0", mode.timestep, frame.tgmm_blob);
    // std::vector<ScaleBlob*> blobs = frame.tgmm_blob->children;
    // // std::vector<ScaleBlob*> blobs = collect_blobs(frame.tgmm_blob, -1, -1);
    // printf("draw %p = %lu blobs.\n", frame.tgmm_blob, blobs.size());
    // // filter.clear();
    // filter.capture(exp->get(mode.timestep));
    // filter.draw_blobs(blobs, ".m");
    // filter.commit(store.buf[0]);
    // return (last_nrrd = store.buf[0]);   
  }
  if(!strcmp(mode.name, "cells")){
    ArFrameData frame = get(mode.timestep);
    std::vector<Cell*> cellsi = cells[int(mode.timestep)];

    filter.clear();
    std::vector<ScaleBlob*> draw;
    std::vector<ScaleBlob*> color;
    for(Cell *cell : cellsi){
      // singleton[0] = cell->blob;
      if(!cell->newcell) draw.push_back(cell->blob);
      else          color.push_back(cell->blob);
    }
    filter.draw_blobs(draw, "lm");
    filter.color_blobs(color, 4.f);
    filter.scale(0.5);
    filter.commit(store.buf[0]);


    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "diff_blobs")){
    ArFrameData frame = get(mode.timestep);

    // float scalemin = frame.scales[mode.blob.scale] - frame.scale_eps;
    // float scalemax = frame.scales[mode.blob.scale] + frame.scale_eps;
    // std::vector<ScaleBlob*> blobs = collect_blobs(frame.blob, scalemin, scalemax);

    std::vector<ScaleBlob*> blobs;
    std::vector<Cell*> cellsi = cells[int(mode.timestep)];
    for(Cell *cell : cellsi){
      blobs.push_back(cell->blob);
    }

    filter.clear();
    filter.draw_blobs(blobs, "gm");
    // filter.scale(0.8);

    filter.difference_image(exp->get(mode.timestep));

    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);   
  }
  if(!strcmp(mode.name, "blobs") || !strcmp(mode.name, "blobs_succs")){
    ArFrameData frame = get(mode.timestep);
    std::vector<Cell*> cellsi = cells[int(mode.timestep)];
    // printf("get scale %d\n", mode.blob.scale);
    // float scalemin = frame.scales[mode.blob.scale] - frame.scale_eps;
    // float scalemax = frame.scales[mode.blob.scale] + frame.scale_eps;
    // printf("mode %s. view scale %d with %.2f %.2f\n", mode.name, scalemin, scalemax);
    // std::vector<ScaleBlob*> blobs = collect_blobs(frame.blob, scalemin, scalemax);
    // std::vector<ScaleBlob*> blobs = collect_blobs(frame.blob, -1, -1);
    
    // draw the blobs in paths at the current timestep.
    // if(!strcmp(mode.geom, "paths")){
    //   for(std::vector<ScaleBlob*> path : mode.highlight.paths){
    //     if(mode.timestep < path.size()){
    //       blobs.push_back(path[mode.timestep]);
    //     }
    //   }
    // }

    // printf("\n");
    filter.clear();
    // for(int i=0;i<mode.highlight.highlight_loci.size();i++){
    //   blobs.push_back(mode.highlight.highlight_loci[i]);
    // }
    // blobs.push_back(mode.highlight.highlight_loci);
    // blobs.insert(blobs.end(), mode.highlight.highlight_loci.begin(), mode.highlight.highlight_loci.end());
    // filter.capture(exp->get(mode.timestep));
    // std::vector<ScaleBlob*> blobs = collect_blobs(frame.blob, scalemin, scalemax);
    // filter.capture(exp->get(mode.timestep));
    // filter.draw_blobs(blobs, "g+");
    std::vector<ScaleBlob*> blobs;
    for(Cell *c : cellsi)blobs.push_back(c->blob);
    filter.draw_blobs(blobs, "g+");
    printf("draw blobs %d\n", mode.highlight.blobs.size());
    filter.draw_blobs(mode.highlight.blobs, ".m");


    // filter.draw_blobs(blobs, ".m");
    // printf("model type = %c\n", blobs[0]->model.type);

    // for(ScaleBlob *sb : mode.highlight.blobs){
    // if(mode.highlight.highlight_loci.size() > 0){
    //   // printf("highlight %lu; scale = %.2f \n", mode.highlight.highlight_loci.size(), mode.highlight.highlight_loci[0]->scale);
    //   // filter.color_blobs(mode.highlight.highlight_loci, 2.f);
    // }
    // if(mode.highlight.highlight_loci.size() > 1){
    //   int i=0;
    //   int j=0;
    //   for(int i=0;i<mode.highlight.highlight_loci.size();++i){
    //     for(int j=0;j<mode.highlight.highlight_loci.size();++j){
    //       printf("  distance %d %d -- %.4f\n", i, j, 
    //         mode.highlight.highlight_loci[i]
    //           ->distance(
    //         mode.highlight.highlight_loci[j])
    //         );
    //     }
    //   }
    // }
    // }

    // show all successors of  blob.
    // if(!strcmp(mode.name, "blobs_succs")){
    //   BSPTree<ScaleBlob> *t1 = &frames[ mode.highlight.timestep - exp->low    ].bspblobs;
    //   std::vector<ScaleBlob*> succs;
    //   for(ScaleBlob *sb : mode.highlight.blobs){
    //     std::vector<ScaleBlob*> potential;
    //     t1->find_within_distance(potential, sb->position, 1000.f);
    //     for(int i=0;i<potential.size();i++){
    //       if(sb->n>1 && potential[i]->n>1 && sb->distance(potential[i]) <= 1.f){
    //         // if(sb->detCov/potential[i]->detCov < 2.f && potential[i]->detCov/sb->detCov < 2.f){
    //           // volume cannot more than double or half.
    //           succs.push_back(potential[i]);
    //           potential[i]->pred.push_back(sb);
    //         // }
    //       }
    //       // printf("%d.", i);
    //     }
    //   }
    //   filter.color_blobs(succs, 4.f);
    // }
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "gaussian")){
    // printf("mode %s.\n", mode.name);
    ArFrameData frame = get(mode.timestep);
    float scale = frame.scales[mode.blob.scale];
    // float scale = 2.f;
    filter.capture(exp->get(mode.timestep));
    DiscreteKernel kernel = filter.gaussian(scale, int(scale*4));
    filter.set_kernel(kernel);
    filter.max1();
    filter.filter();
    // filter.highlight(filter.find_maxima());
    filter.commit(store.buf[0]);
    kernel.destroy();
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "laplacian")){
    
    ArFrameData frame = get(mode.timestep);

    float scale = frame.scales[mode.blob.scale];
    filter.capture(exp->get(mode.timestep));
    DiscreteKernel kernel1 = filter.gaussian(scale, int(scale*4));

    filter.max1();    
    filter.set_kernel(kernel1);
    filter.filter();
    filter.laplacian3d();
    filter.normalize();

    kernel1.destroy();
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "hessian")){
    
    ArFrameData frame = get(mode.timestep);

    float scale = frame.scales[mode.blob.scale];
    filter.capture(exp->get(mode.timestep));
    DiscreteKernel kernel1 = filter.gaussian(scale, int(scale*4));

    filter.max1();    
    filter.set_kernel(kernel1);
    filter.filter();
    filter.hessian3d();
    filter.normalize();

    kernel1.destroy();
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "masked")){
    filter.capture(exp->get(mode.timestep));
    ArFrameData frame = get(mode.timestep);
    float scale = frame.scales[mode.blob.scale];
    filter.max1();
    filter.lapofgaussian_masked(scale, false);
    filter.normalize();
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "maxima")){
    
    ArFrameData frame = get(mode.timestep);

    float scale = frame.scales[mode.blob.scale];
    filter.capture(exp->get(mode.timestep));
    DiscreteKernel kernel1 = filter.gaussian(scale, int(scale*4));

    filter.max1();    
    filter.set_kernel(kernel1);
    filter.filter();
    filter.laplacian3d();
    filter.normalize();

    filter.highlight(filter.find_maxima());

    kernel1.destroy();
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "levelsets")){
    Nrrd* nrrd = exp->get(mode.timestep);
    float *data = (float*)nrrd->data;
    filter.capture(nrrd);
    ArFrameData frame = get(mode.timestep);
    float scale = frame.scales[mode.blob.scale];
    filter.max1();
    // filter.lapofgaussian_masked(scale, 0);

    // filter.max1();

    DiscreteKernel kernel = filter.gaussian(scale, int(scale*4));
    filter.set_kernel(kernel);
    kernel.destroy();
    filter.filter();
    filter.laplacian3d();

    filter.normalize();
    filter.show_blobs(1);
    // filter.highlight(filter.find_maxima());
    // filter.normalize(3.f);
    filter.commit(store.buf[0]);
    for(int i=0;i<filter.self.w3;i++){
      if(((float*)(store.buf[0]->data))[i] > 0){
        // float *out = 
        ((float*)(store.buf[0]->data))[i] = data[i];
      }
    }
    return (last_nrrd = store.buf[0]);
  }
  if(!strcmp(mode.name, "show_blobs")){
    filter.capture(exp->get(mode.timestep));
    ArFrameData frame = get(mode.timestep);
    float scale = frame.scales[mode.blob.scale];
    filter.max1();
    filter.lapofgaussian_masked(scale, 0);
    filter.normalize();
    filter.show_blobs(0);
    filter.normalize();
    filter.commit(store.buf[0]);
    return (last_nrrd = store.buf[0]);
  }
}
/*
 * emit the current representation
 * as a set of .nrrd files.         ASSUMES 3 DIMENSIONS
 */
void ArPipeline::emit(ReprMode &mode, std::string suffix, int low, int high){
  // printf("emit %d %d\n", low, high);
  for(int i=low; i<=high; i++){
    if(i < exp->low)continue;
    if(i > exp->high)continue;
    mode.timestep = i;
    printf("emit to %s\n", (exp->getfilepath(i)+suffix).c_str());
    Nrrd *nrrdfloat = repr(mode);
    Nrrd *nrrdshort = nrrdNew();
    nrrdConvert(nrrdshort, nrrdfloat, nrrdTypeUShort);
    int size = nrrdshort->axis[0].size * nrrdshort->axis[1].size * nrrdshort->axis[2].size;
    short *datashort = (short*)nrrdshort->data;
    float *datafloat = (float*)nrrdfloat->data;
    for(int i=0;i<size;++i){
      datashort[i] = (unsigned short) (datafloat[i] * 32767.0);
    }
    // nrrdCopy(out, asfloat);
    nrrdSave((exp->getfilepath(i)+suffix).c_str(), nrrdshort, NULL);
    nrrdNuke(nrrdshort);
  }
}

#define WRITET(T, x) {T xx = (T)x; fwrite(&(xx), sizeof(T), 1, file);}
#define WRITE(x)    {fwrite(&(x), sizeof(x), 1, file);}
#define READ(x) (good&=!!(fread(&(x), sizeof(x), 1, file)));

void ArPipeline::loadframe(int timestep){
  int good = 1;
  std::string path = (exp->getfilepath(timestep)+".blobs");
  FILE *file = fopen(path.c_str(),"rb");
  if(!file)return;
  // printf("load %s\n", path.c_str());
  int nblobs = 1;
  
  std::map<void*, ScaleBlob*> label_to_blob;
  label_to_blob[0] = 0;

  ScaleBlob *root;
  std::vector<ScaleBlob*> blobs;
  while(nblobs){
    // printf("n=%d\n", nblobs);
    ScaleBlob *blob = new ScaleBlob;
    blobs.push_back(blob);
    void *label;
    int nchildren;
    READ(label);
    label_to_blob[label] = blob;
    READ(nchildren);
    // printf("n += %d\n", nchildren);
    nblobs += nchildren;
    for(int i=0;i<nchildren;i++){
      ScaleBlob *child;
      READ(child);
      blob->children.push_back(child);
    }
    READ(blob->parent);
    READ(blob->mode);
    READ(blob->imode);
    READ(blob->position);
    READ(blob->shape);
    READ(blob->timestep);
    READ(blob->covariance);
    READ(blob->invCov);
    READ(blob->detCov);
    READ(blob->pdfCoef);
    READ(blob->min);
    READ(blob->max);
    READ(blob->scale);
    READ(blob->n);
    READ(blob->npass);
    READ(blob->peakvalue);
    READ(blob->GMM.nk);
    READ(blob->GMM.alpha);
    READ(blob->GMM.gaus_scale);
    nblobs -= 1;
  }
  fclose(file);

  if(!good)return;

  // convert labels into real pointers
  for(ScaleBlob *blob : blobs){
    for(int i=0;i<blob->children.size();i++){
      blob->children[i] = label_to_blob[blob->children[i]];
    }
    blob->parent = label_to_blob[blob->parent];
  }

  root = blobs[0];
  ArFrameData &frame = get(timestep);
  frame.complete = true;
  frame.blob = root;
  frame.bspblobs = filter.get_bsp(10);
  for(ScaleBlob* blob : blobs)frame.bspblobs.insert(blob, blob->mode);
  frame.scales = collect_scales(root);
  frame.scale_eps = compute_epsilon(frame.scales);
}
void ArPipeline::saveframe(int timestep){
  std::string path = (exp->getfilepath(timestep)+".blobs");
  printf("write to %s\n", path.c_str());
  FILE *file = fopen(path.c_str(),"wb");
  if(!file){
    printf("cannot write to file.");
  }
  std::queue<ScaleBlob*> blobs;
  blobs.push(get(timestep).blob);
  if(!get(timestep).blob){
    printf("frames(%d).blob = null\n", timestep);
    return;
  }
  while(!blobs.empty()){
    ScaleBlob *sb = blobs.front();
    blobs.pop();
    WRITET(ScaleBlob*, sb);
    WRITET(int, sb->children.size());
    for(ScaleBlob *s : sb->children){
      WRITE(s);
      blobs.push(s);
    }
    WRITE(sb->parent);
    WRITE(sb->mode);
    WRITE(sb->imode);
    WRITE(sb->position);
    WRITE(sb->shape);
    WRITE(sb->timestep);
    WRITE(sb->covariance);
    WRITE(sb->invCov);
    WRITE(sb->detCov);
    WRITE(sb->pdfCoef);
    WRITE(sb->min);
    WRITE(sb->max);
    WRITE(sb->scale);
    WRITE(sb->n);
    WRITE(sb->npass);
    WRITE(sb->peakvalue);
    WRITE(sb->GMM.nk);
    WRITE(sb->GMM.alpha);
    WRITE(sb->GMM.gaus_scale);
  }
  fclose(file);
}


  /* Save/Load processed pipeline. format:
   * filename as [filepath_0].pipeline
   * [int] number of frames processed sequentially from 0.
   * <list of blobs:
   *    [int] address of blob
   *    ... data ...
   * 
   * >
   */

void ArPipeline::save(){


  std::string path0 = exp->getfilepath(exp->low);
  std::replace(path0.begin(), path0.end(), '/', '-');
  path0 = "../rsc/store/s" + path0;

  printf("writing to : %s\n", path0.c_str());

  FILE *file = fopen((path0 + ".pipeline").c_str(),"wb");

  // fwrite("hello!",sizeof(char),6,file);

  // count number of frames that are processed.
  int nframes = 0;
  for(int i=0;i<frames.size();i++){
    if(!frames[i].complete){
      break;
    }else{
      nframes = i+1;
    }
  }

  WRITET(int, nframes);
  // printf("write %d processed\n", nframes);

  for(int i=0;i<nframes;++i){
    // fwrite("ff", sizeof(char), 2, file);
    // write address of root scaleblob.
    WRITET(ScaleBlob*, frames[i].blob);
    // printf("write %p\n", frames[i].blob);
    // fwrite(&frames[i].blob, sizeof(ScaleBlob*), 1, file);
    std::vector<ScaleBlob*> allblobs = frames[i].bspblobs.as_vector();
    
    // write total number of blobs in this frame.
    WRITET(int, allblobs.size());
    // fwrite(&nblobs, sizeof(int), 1, file);

    for(ScaleBlob *sb : allblobs){
      // write each blob
      WRITET(ScaleBlob*, sb);

      WRITET(int, sb->children.size());
      for(ScaleBlob *s : sb->children){
        WRITE(s);
      }

      WRITET(int, sb->pred.size());
      for(ScaleBlob *s : sb->pred){
        WRITE(s);
      }

      WRITET(int, sb->succ.size());
      for(ScaleBlob *s : sb->succ){
        WRITE(s);
      }
      WRITE(sb->parent);
      WRITE(sb->mode);
      WRITE(sb->imode);
      WRITE(sb->position);
      WRITE(sb->shape);
      WRITE(sb->timestep);
      WRITE(sb->covariance);
      WRITE(sb->invCov);
      WRITE(sb->detCov);
      WRITE(sb->pdfCoef);
      WRITE(sb->min);
      WRITE(sb->max);
      WRITE(sb->scale);
      WRITE(sb->n);
      WRITE(sb->npass);
      WRITE(sb->peakvalue);
      // model:
      WRITE(sb->model.alpha);
      WRITE(sb->model.beta);
      WRITE(sb->model.kappa);
      WRITE(sb->model.type);
      WRITE(sb->model.n);
      WRITE(sb->model.error);
      WRITE(sb->model.min);
      WRITE(sb->model.max);
    }
  }

  fclose(file);

  file = fopen((path0 + ".paths").c_str(),"wb");
  printf("write %lu paths\n", paths.size());
  WRITET(int, paths.size());
  for(std::vector<ScaleBlob*> path : paths){
    WRITET(int, path.size());
    for(ScaleBlob* s : path){
      WRITE(s);
    }
  }

  fclose(file);

  file = fopen((path0 + ".paths.txt").c_str(), "w");
  for(std::vector<ScaleBlob*> path : paths){
    for(ScaleBlob* blob : path){
      // if(blob)
      fprintf(file, "%.2f %.2f %.2f; ", blob->position.x, blob->position.y, blob->position.z);
    }
    fprintf(file, "\n");
  }

}
void ArPipeline::load(){
//   printf("loading...\n");
  std::string path0 = exp->getfilepath(exp->low);
  std::replace(path0.begin(), path0.end(), '/', '-');
  path0 = "../rsc/store/s" + path0;

  std::string path0pipeline = path0 + ".pipeline";
  std::string path0paths    = path0 + ".paths";

  if(access(path0pipeline.c_str(), F_OK) == -1)return;
  FILE *file = fopen(path0pipeline.c_str(),"rb");

  char buf[2];
  int good=1;
  int nframes;
  READ(nframes);
  printf("nframes = %d\n", nframes);
  std::vector<ScaleBlob*> frameroots;
  std::unordered_map<ScaleBlob*, ScaleBlob*> allblobs;
  std::vector<std::vector<ScaleBlob*>> frameblobs;
  for(int i=0;i<nframes;++i){
    ScaleBlob *rootaddr;
    int   nblobs;
    READ(rootaddr);
    READ(nblobs);
    frameroots.push_back(rootaddr);
    // printf("root %p has %d blobs.\n", rootaddr, nblobs);
    std::vector<ScaleBlob*> blobs;
    for(int i=0;i<nblobs;++i){
      ScaleBlob* p0;
      ScaleBlob* label;
      int nchildren, npred, nsucc;
      READ(label);
      if(!label)continue;
      ScaleBlob *blob = new ScaleBlob();
      READ(nchildren);
      for(int i=0;i<nchildren;++i){
        READ(p0);
        blob->children.push_back((ScaleBlob*)p0);
      }
      READ(npred);
      for(int i=0;i<npred;++i){
        READ(p0);
        blob->pred.push_back((ScaleBlob*)p0);
      }
      READ(nsucc);
      for(int i=0;i<nsucc;++i){
        READ(p0);
        blob->succ.push_back((ScaleBlob*)p0);
      }
      READ(blob->parent);
      READ(blob->mode);
      READ(blob->imode);
      READ(blob->position);
      READ(blob->shape);
      READ(blob->timestep);
      READ(blob->covariance);
      READ(blob->invCov);
      READ(blob->detCov);
      READ(blob->pdfCoef);
      READ(blob->min);
      READ(blob->max);
      READ(blob->scale);
      READ(blob->n);
      READ(blob->npass);
      READ(blob->peakvalue);
      // model:
      READ(blob->model.alpha);
      READ(blob->model.beta);
      READ(blob->model.kappa);
      READ(blob->model.type);
      READ(blob->model.n);
      READ(blob->model.error);
      READ(blob->model.min);
      READ(blob->model.max);

      // printf("read blob->n = %.2f\n", blob->n);

      allblobs[label] = blob;
      // allblobs.insert({label, blob});
      // printf("%d: %p\t%p\n", allblobs.size(), label, blob);
      blobs.push_back(blob);
    }
    frameblobs.push_back(blobs);
  }
  // allblobs[0] = 0;
  if(!good){
    fprintf(stderr, "pipeline.load: read error.\n");
    return;
    // exit(0);
  }
  // int itr=0;
  // for(auto  it = allblobs.begin(); it != allblobs.end(); ++it){
  //   printf("%d %p %p\n", ++itr, it->first, it->second);
  // }
  // printf("fixing...\n");
  int itr=0;
  // printf("size=%d\n", allblobs.size());
  for(auto  elt = allblobs.begin(); elt != allblobs.end(); ++elt){
    if(!elt->second)continue;
    // printf("%d %p %p\n", ++itr, elt->first, elt->second);
    // printf("parents...");
    // printf("%p; \n", elt->second->parent);
    if(elt->second->parent)
      elt->second->parent = allblobs[elt->second->parent];
    // printf("children...");
    for(int i=0;i<elt->second->children.size(); ++i){
      // printf("%d/%d; ", i, elt->second->children.size());
      elt->second->children[i] = allblobs[elt->second->children[i]];
    }
    // printf("preds...");
    for(int i=0;i<elt->second->pred.size(); ++i){
      // printf("%d/%d; ", i, elt->second->pred.size());
      elt->second->pred[i] = allblobs[elt->second->pred[i]];
    }
    // printf("succs...");
    for(int i=0;i<elt->second->succ.size(); ++i){
      // printf("%d/%d; ", i, elt->second->succ.size());
      elt->second->succ[i] = allblobs[elt->second->succ[i]];
    }
  }
  printf("\nread......\n");
  // printf("all blobs:");
  int i=0;
  for(std::pair<ScaleBlob*, ScaleBlob*> elt : allblobs){
    // printf("%d: %p %p\n", i, elt.first, elt.second);
    ++i;
  }
  // printf("pushing\n");
  // printf("frames %d %d\n", nframes, frames.size());
  for(int i=0;i<nframes;i++){
    ArFrameData fd;
    fd.blob      = allblobs[frameroots[i]];
    // printf("%p\n", fd.blob);
    fd.scales    = collect_scales(fd.blob);
    fd.scale_eps = compute_epsilon(fd.scales);
    fd.complete  = true;
    fd.bspblobs = filter.get_bsp(10);
    
    BSPTree<ScaleBlob> *bsptree = &fd.bspblobs;
    for(ScaleBlob *sb : frameblobs[i]){
      bsptree->insert(sb, sb->position);
    }

    // DEBUG

    // for(float f : fd.scales){
    //   printf("scale %.2f;\n", f);
    // }
    // /////
    if(i == frames.size())frames.push_back(fd);
    else frames[i] = fd;
  }

  // load paths...

  if(access(path0paths.c_str(), F_OK) == -1)return;
  file = fopen(path0paths.c_str(),"rb");

  int n;
  READ(n);
  printf("read %d paths.\n", n);
  for(int i=0;i<n;++i){
    std::vector<ScaleBlob*> path;
    int pathlen;
    ScaleBlob *blob;
    READ(pathlen);
    for(int i=0;i<pathlen;++i){
      READ(blob);
      path.push_back(allblobs[blob]);
    }
    paths.push_back(path);
  }

  // loadTGMM();

  // done.
}

#undef WRITE
#undef WRITET
#undef READ

// void ArPipeline::getTGMMOutput(std::string path, int low, int high){
//   printf("load tgmm output\n");
// }
static void print_element_names(xmlNode * a_node)
 {
    xmlNode *cur_node = NULL;

    for (cur_node = a_node; cur_node; cur_node =
         cur_node->next) {
       if (cur_node->type == XML_ELEMENT_NODE) {
          printf("node type: Element, name: %s\n",
               cur_node->name);
       }
       print_element_names(cur_node->children);
    }
 }

static vec3 from_string(char *in){
  char *token;
  vec3 ret(0,0,0);
 
  token = strtok(in, " ");
  ret.x = std::stof(token);
  token = strtok(NULL, " ");
  ret.y = std::stof(token);
  token = strtok(NULL, " ");
  ret.z = std::stof(token);
  // printf("token = %s", token);
  return ret;
}
static dmat3x3 dmat3_from_string(char *W){
  return mat3x3(5,0,0, 0,5,0, 0,0,5);
}
void ArPipeline::loadTGMM(){
  xmlDoc *doc = NULL;
  xmlNode *root = NULL;


  // map of blobs between frames
  std::map<int, ScaleBlob*> maps[] = {
    std::map<int, ScaleBlob*>(),
    std::map<int, ScaleBlob*>()
  };int mapi = 0;

  for(int i=exp->low, ii=0; i<=exp->high;i++, ii++){
    // printf("TGMM load %s.\n", exp->gettgmmpath(i).c_str());
    // doc = xmlReadFile(exp->gettgmmpath(i).c_str(), NULL, 0);
    // if(!doc){
    //   printf("could not read file.\n");
    //   continue;
    // }else{
    //   printf("read %s\n", exp->gettgmmpath(i).c_str());
    // }

    root = xmlDocGetRootElement(doc);


    xmlNode *cur = 0;
    
    // frames[ii].tgmm_blob = new ScaleBlob();
    maps[mapi].clear();

    int framen = 0;

    for(cur = root->children; cur; cur = cur->next){
      if (cur->type == XML_ELEMENT_NODE){
        // GaussianMixtureModel
        ScaleBlob *blob = new ScaleBlob();
        xmlChar *m = xmlGetProp(cur, (xmlChar*)"m");
        xmlChar *W = xmlGetProp(cur, (xmlChar*)"W");
        xmlChar *s_id = xmlGetProp(cur, (xmlChar*)"id");
        xmlChar *s_parent = xmlGetProp(cur, (xmlChar*)"parent");

        int id = std::stof((char*) s_id);
        int pred = std::stof((char*) s_parent);

        maps[mapi][id] = blob;
        if(pred != -1){
          maps[mapi][id]->pred.push_back(maps[!mapi][pred]);
          maps[!mapi][pred]->succ.push_back(maps[mapi][id]);
        }
        // printf("m=%s\n", m);
        // printf("W=%s\n", W);
        blob->n = 100;
        blob->scale = 0;
        blob->position = from_string((char*)m);
        blob->min = blob->position - dvec3(10,10,10);
        blob->max = blob->position + dvec3(10,10,10);
        blob->shape = dmat3_from_string((char*)W);
        // printf("m = %.2f %.2f %.2f\n", blob->position.x, blob->position.y, blob->position.z);
        // frames[ii].tgmm_blob->children.push_back(blob);
        // printf("nchildren = %d\n", frames[ii].tgmm_blob->children.size());
      }
    }
    mapi = !mapi;
    // printf("size[%d] = %p = %d\n", ii, frames[ii].tgmm_blob, frames[ii].tgmm_blob->children.size());
    // print_element_names(root);
    xmlFreeDoc(doc);       // free document
    xmlCleanupParser();    // Free globals
  }

}