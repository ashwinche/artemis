#ifndef UTIL_H
#define UTIL_H

#include <chrono>

static void tick(std::string text){
  using namespace std;
  using namespace std::chrono; 
  static auto clock = high_resolution_clock::now();
  auto next = high_resolution_clock::now();
  long d = duration_cast<milliseconds>(next-clock).count();
  printf("%s",text.c_str());
  printf(">> elapsed %lums\n\n",d);
  clock = next;
}
#endif