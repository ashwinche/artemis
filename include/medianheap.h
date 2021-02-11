#ifndef MEDIAN_HEAP
#define MEDIAN_HEAP

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <cstdlib>

#define LCHILD(x) ((x<<1) + 1)
#define RCHILD(x) ((x<<1) + 2)
#define PARENT(x) ((x - 1) / 2)
#define SWAPELTS(t1, t2, heap, x, y) (t1=heap[x],t2=index[x],heap[x]=heap[y],index[x]=index[y],heap[y]=t1,index[y]=t2)

#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)>(y))?(x):(y))


typedef struct{
  int value;
  int index;
}Node;
class MedianHeap{
private:
  Node *heaps;
public:
  MedianHeap(int n){
    this->n = n;
    minl = n/2;
    maxl = n-minl;  // maxl might be bigger than minl.

    heaps = new Node[n];
    memset(heaps,0,n*sizeof(int));
    minh = heaps;
    maxh = heaps+minl;

    index = new int[n];
    for(int i=0;i<n; ++i){
      heaps[i].index=i;
      index[i]=i;
    }
  }
  int n;
  int  *index;
  Node *minh, *maxh;
  int   minl,  maxl;
  int median(){
    return maxh[0].value;
  }
  void assert(bool b, char *s, char *name){
    if(!b){
      printf("failed %s: %s\n", s, name);
      print();
      exit(1);
    }
  }
  void verify(char *s){
    for(int i=0;i<minl;++i){
      assert(minh[i].value >= maxh[0].value, s, "minheap is greater than median.");
      if(i){
        assert(minh[i].value >= minh[PARENT(i)].value, s, "minheap property");
      }
    }
    for(int i=0;i<maxl;++i){
      assert(maxh[i].value <= maxh[0].value, s, "maxheap is less than median.");
      if(i){
        assert(maxh[i].value <= maxh[PARENT(i)].value, s, "maxheap property");
      }
    }
  }
  void printall(){
    for(int i=0;i<n;++i){
      printf("%4d ",heaps[i]);
    }printf("\n");
  }
  void print(){
    printf("\n\n\t");
    // printf("~~~\n");
    for(int i=0;i<minl;++i){
      printf("%d\t",minh[i].index);
    }
    printf("\n\t");
    for(int i=0;i<minl;++i){
      printf("%d\t",minh[i].value);
    }
    printf("\n\n\t\t");
    printf(" > %d > ", median());
    printf("\n\n\t");
    for(int i=0;i<maxl;++i){
      printf("%d\t",maxh[i].index);
    }
    printf("\n\t");
    for(int i=0;i<maxl;++i){
      printf("%d\t",maxh[i].value);
    }
    // printf("\n~~~\n");
    printf("\n\n");
    printf("\n");
  }
  inline int swap_elts(Node *const heap, const int xi, const int yi){
    Node t1  = heap[xi];
    heap[xi] = heap[yi];
    heap[yi] = t1;

    int t2 = index[heap[xi].index];
    index[heap[xi].index] = index[heap[yi].index];
    index[heap[yi].index] = t2;
  }
  inline int swap(int ix, int iy){
    Node t = heaps[ix];
    heaps[ix] = heaps[iy];
    heaps[iy] = t;

    int ti = index[heaps[ix].index];
    index[heaps[ix].index] = index[heaps[iy].index];
    index[heaps[iy].index] = ti;
  }
  void replace(int ind, int value){

    // printf("%3d: %3d -> %3d\n", ind, heaps[index[ind]], value);
    // verify();
    Node *p;
    int pi, xi;
    int pv, xv;
    int t1, t2;

    Node x{value, ind};

    bool fixed = false;

    xi = index[ind];
    xv = value;

    static Node MH_MAX{INT_MAX, 0};
    static Node MH_MIN{INT_MIN, 0};

    // printf("insert (%d,%d). index: %d\n",ind,value,xi);
    // verify("0");
    // bubble up.

    if(xi < minl){  // minheap.
      // printf("minheap\n");
      while(xi){
        // move xi upward until it hits the top, if needed.
        pi = PARENT(xi);
        pv = minh[pi].value;
        if(xv < pv){
          minh[xi] = minh[pi];
          index[minh[xi].index] = xi;
          fixed = true;
        }else{
          minh[xi] = x;
          index[minh[xi].index] = xi;
          break;
        }
        xi = pi;
      }if(!xi){
        minh[xi] = x;
      }
      if(!xi && x.value < maxh[0].value){
        // printf("before swap\n");
        // print();
        swap(0, minl);  // move x to the maxheap.
        xi = minl;      // signal to percolate down in the maxheap.
      }else if(fixed){
        // verify("2");
        return;
      }
    }else{          // maxheap.
      xi = xi-minl;
      // printf("maxheap\n");
      // printf("xi=%d\n",xi);
      // print();
      while(xi){
        // move xi upward until it hits the top, if needed.
        pi = PARENT(xi);
        pv = maxh[pi].value;
        if(xv > pv){
          maxh[xi] = maxh[pi];
          index[maxh[xi].index] = minl + xi;
          fixed = true;
        }else{
          maxh[xi] = x;
          index[maxh[xi].index] = minl + xi;
          break;
        }
        xi = pi;
      }if(!xi){
        maxh[xi] = x;
      }
      if(!xi && x.value > minh[0].value){
        // printf("swap\n");
        // print();
        swap(0, minl);  // move x to the minheap.
        xi = 0;         // signal to percolate down in the minheap.
      }else if(fixed){
        // verify("4");
        return;
      }else xi = xi+minl;
    }

    // print();
    // printf("bubble down\n");

    // bubble down.

    int li, ri;
    Node *l, *r;

    bool xl, xr, lr;

    if(xi < minl){  // minheap
      // printf("bubble down minheap\n");
      // print();
      // check children
      for(;;){
        li=LCHILD(xi);
        ri=RCHILD(xi);

        // prevent overflow when there are
        // no more children.
        l=(li<minl)?(minh+li):&MH_MAX;
        r=(ri<minl)?(minh+ri):&MH_MAX;
        
        // used to correctly swap elements.
        xl=(xv       <= l->value);
        xr=(xv       <= r->value);
        lr=(l->value <= r->value);
        
        if(xl && xr){
          minh[xi] = x;
          index[minh[xi].index] = xi;
          // verify("6");
          return;
        }else if(!xl && lr){
          // swap with left.
          minh[xi] = minh[li];
          index[minh[xi].index] = xi;
          xi = li;
        }else{
          // swap with right.
          minh[xi] = minh[ri];
          index[minh[xi].index] = xi;
          xi = ri;
        }
      }
    }else{
      xi -= minl;
      // printf("bubble down maxheap\n");
      // print();
      for(;;){
        li=LCHILD(xi);
        ri=RCHILD(xi);

        // prevent overflow when there are
        // no more children.
        l=(li<maxl)?(maxh+li):&MH_MIN;
        r=(ri<maxl)?(maxh+ri):&MH_MIN;
        
        // used to correctly swap elements.
        xl=(xv       <= l->value);
        xr=(xv       <= r->value);
        lr=(l->value <= r->value);
        
        if(!xl && !xr){
          // printf(" stay\n");
          maxh[xi] = x;
          index[maxh[xi].index] = minl + xi;
          // verify("7");
          return;
        }else if(xl && !lr){
          // printf("left ");
          // swap with left.
          maxh[xi] = maxh[li];
          index[maxh[xi].index] = minl + xi;
          xi = li;
        }else{
          // printf("right ");
          // swap with right.
          maxh[xi] = maxh[ri];
          index[maxh[xi].index] = minl + xi;
          xi = ri;
        }
      }
    }
    // verify("8");


  }
};

#undef LCHILD
#undef RCHILD
#undef PARENT

#undef min
#undef max
#endif