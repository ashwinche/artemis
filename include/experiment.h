#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <teem/nrrd.h>
#include <string>
#include <map>
#include "filter.h"

struct NrrdFrame{
  float n;
  long accessed;
  std::string path;
  Nrrd* nrrd;
};

// Exposes methods for interacting with
// the raw data of the experiment.
class ArExperiment{
public:
  ArExperiment(std::string path, int low, int high, int mem_cap);
  NrrdFrame* frames;
  Nrrd* get(float n, bool force = false);
  Nrrd* copy(int n);
  void interpolate(Nrrd* f0, Nrrd *f1, Nrrd *fx, float alpha);
  std::string getfilepath(int n);
  // std::string gettgmmpath(int n);
  int low;
  int high;
  std::string filepath;
private:
  ArFilter filter;
  std::string *paths;
  // std::string *tgmmpaths;
  int nframes;
  int npaths;
  long time;
};

#endif