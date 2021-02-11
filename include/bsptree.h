#ifndef BSPTREE_H
#define BSPTREE_H

#include "blob.h"
#include <glm/glm.hpp>
#include <vector>
#include <algorithm>

using namespace glm;

struct plane{
  vec3 point;
  vec3 normal;
};
template<class T>
class BSPTree{
private:
  struct E{
    T *v;
    vec3 p;
  };
public:
  struct D{
    D(T *v, float d) : v(v), d(d){}
    T *v;
    float d;
  };
  BSPTree(BSPTree *parent, int depth, vec3 min, vec3 max) : parent(parent), min(min), max(max), children(0), n(0){
    split(depth);
  }
  BSPTree() : min(vec3(0)), max(vec3(0)), parent(0), children(0), n(0){}
  
  vec3 min;
  vec3 max;
  plane cutting;
  BSPTree *parent;
  BSPTree *children;      // enforce: array of size either 0 or 2.
  std::vector<E> members;
  int n;                  // number of children = 

  void indent(int n){
    for(int i=0;i<n;i++)printf(" ");
  }
  void print(int c = 0){
    indent(c);
    printf("T%d (%.1f %.1f %.1f)  %.0f %.0f %.0f  (%.1f %.1f %.1f)\n",c, min.x,min.y,min.z, cutting.normal.x, cutting.normal.y,cutting.normal.z, max.x, max.y, max.z);
    if(children){
      children[0].print(c+1);
      children[1].print(c+1);
    }else{
      for(int i=0;i<members.size();i++){
        indent(c+1);
        printf("(%.1f %.1f %.1f) x\n", members[i].p.x, members[i].p.y, members[i].p.z);
      }
    }
  }
  void find_within_distance(std::vector<T*> &list, vec3 q, float maxdsq){
    // check should ignore box? too far from q.
    if(distsq_point_bbox(q) > maxdsq){
      return; 
    }
    for(E e : members){
      float edsq = distsq(q, e.p);
      if(edsq < maxdsq){
        list.push_back(e.v);
      }
    }
    if(children){
      children[0].find_within_distance(list, q, maxdsq);
      children[1].find_within_distance(list, q, maxdsq);
    }
  }
  T *find_closest(vec3 q){
    T* ret = 0;
    float dist = 0;
    find_closest(q, ret, dist);
    return ret;
  }
  void find_closest(vec3 q, T* &v, float &maxdsq){
    // printf("find closest to %.0f %.0f %.0f in [%.1f %.1f] [%.1f %.1f] [%.1f %.1f].\n", q.x, q.y, q.z, min.x, max.x, min.y, max.y, min.z, max.z);
    if(!v){
      // printf("v=%p\n",v);
      // no candidate point found yet. search for one.
      if(children){
        children[0].find_closest(q, v, maxdsq);
        children[1].find_closest(q, v, maxdsq);
      }else{
        for(E e : members){
          // printf("leaf.\n");
          float dsq = distsq(q, e.p);
          if(!v || dsq < maxdsq){
            maxdsq = dsq;
            v = e.v;
          }
        }
      }
    }
    // printf("v=%p!\n",v);
    // if q is within the appropriate distance to this cell...
    if(distsq_point_bbox(q) < maxdsq){
      if(children){
        children[0].find_closest(q, v, maxdsq);
        children[1].find_closest(q, v, maxdsq);
      }else{
        for(E e : members){
          float dsq = distsq(q, e.p);
          if(dsq < maxdsq){
            maxdsq = dsq;
            v = e.v;
          }
        }
      }
    }
  }
  // distance between p and this BSPTree's bounding box (given by min/max)
  float distsq_point_bbox(vec3 p){
    vec3 q(p);  // closest point within box to p.
    
    // constrain q to be within the box.
         if(q.x<min.x)q.x = min.x;
    else if(q.x>max.x)q.x = max.x;
    
         if(q.y<min.y)q.y = min.y;
    else if(q.y>max.y)q.y = max.y;
    
         if(q.z<min.z)q.z = min.z;
    else if(q.z>max.z)q.z = max.z;
    
    return distsq(p,q);
  }
  std::vector<T*> as_vector(){
    std::vector<T*> ret;
    for(int i=0;i<members.size();i++){
      // printf("i=%d\n", members.size());
      ret.push_back(members[i].v);
    }
    if(children){
      std::vector<T*> lc = children[0].as_vector();
      std::vector<T*> rc = children[1].as_vector();
      ret.insert(ret.end(), lc.begin(), lc.end());
      ret.insert(ret.end(), rc.begin(), rc.end());
    }
    return ret;
  }

  void clear(){
    members.clear();
    if(children){
      children[0].clear();
      children[1].clear();
    }
  }
  void split(int depth){
    if(depth == 0)return;

    children = new BSPTree<T>[2];

    float dx = max.x - min.x;
    float dy = max.y - min.y;
    float dz = max.z - min.z;

    vec3 center = 0.5f*(min+max);

    cutting.point = center; // splitting plane intersects center.
    
    if(dx >= dy && dx >= dz)      cutting.normal = vec3(1,0,0);
    else if(dx >= dy && dx >= dz) cutting.normal = vec3(0,1,0);
    else cutting.normal = vec3(0,0,1);

    children[0] = BSPTree(this, depth-1, min, max - cutting.normal*((max-min)/2.f));
    children[1] = BSPTree(this, depth-1, min + cutting.normal*((max-min)/2.f), max);

    for(int i=0;i<members.size();i++){
      insert(members[i].v, members[i].p);
    }
    members.clear();
  }
  static float distsq(vec3 a, vec3 b){
    vec3 v = a-b;
    return v.x*v.x + v.y*v.y + v.z*v.z;
  }

  BSPTree *getsubtree(vec3 p){
    if(dot(p - cutting.point,cutting.normal) < 0){
      return children;
    }else{
      return children+1;
    }
  }
  void insert(T *v, vec3 p){
    n += 1;
    if(!children){
      E elt;
      elt.v = v;
      elt.p = p;
      members.push_back(elt);
      return;
    }
    getsubtree(p)->insert(v,p);
  }
  bool remove(T *v){
    for(int i=0;i<members.size();i++){
      if(v[i] == v){
        members.erase(members.begin()+i);
        return true;
      }
    }
    if(children && children[0]->remove(v))return true;
    if(children && children[1]->remove(v))return true;
    return false;
  }
  bool contains(T *v){
    if(std::find(members.begin(), members.end(),v) != members.end())return true;
    if(children && children[0]->contains(v))return true;
    if(children && children[1]->contains(v))return true;
    return false;
  }
};

#endif