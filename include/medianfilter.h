
#include <vector>
#include <stdio.h>
#include <cstdlib>

// #include "red_black_tree.h"
#include "medianheap.h"
// #include "quickselect.h"

template<int R>
void median_filter_3D(short *data, int *dims, short *out, int *fmin, int *fsiz, int *fmax){
  int verbose = 0;
  setbuf(stdout, NULL);

  // ensure fmin + fsiz = fmax;
  for(int i=0;i<3;++i){
    if(fmin[i] + fsiz[i] != fmax[i]){
      fprintf(stderr,"malformed argument: fmin[%d] + fsiz[%d] =/= fmax[%d]\n",i,i,i);
      exit(0);
    }
  }
  
  // implement Perreault's algorithm using trees.
  const int N = R*2+1;
  // RBRTree *frame = (RBRTree*)malloc(fsiz[0]*fsiz[1]*sizeof(RBRTree));
  MedianHeap window(N*N*N);

  printf("init window\n");
  // 1. initialize window at origin corner.
  int ind = 0;
  for(int z=fmin[2];z<fmin[2]+N;++z){
    for(int y=fmin[1];y<fmin[1]+N;++y){
      for(int x=fmin[0];x<fmin[0]+N;++x){
        int v = data[z*dims[0]*dims[1] + y*dims[0] + x];
        window.replace(ind++,v);
      }
    }
  }
  verbose = 0;
    
  // current point to be medianed
  int x = fmin[0] + R;
  int y = fmin[1] + R;
  int z = fmin[2] + R;
  
  ind = z*dims[0]*dims[1] + y*dims[0] + x;

  // direction of window's travel, x and y.
  int dir_x = 1;
  int dir_y = 1;

  // offset with which data is stored in window
  int off_x = 0;
  int off_y = 0;
  int off_z = 0;

  // position of window within frmae
  int win_x = 0;
  int win_y = 0;

  // depth of window relative to start
  int depth_z = 0;

  int total = (dims[0]-2*R)*(dims[1]-2*R)*(dims[2]-2*R);
  int iteration = 0;


  for(;;){

    // printf("(%d %d %d %d);",x,y,z,ind);
    // for(int j=0;j<N*N;++j){
    //   window[j]->printall();
    // }
    verbose = 0;
    // dprint("index %d...\n", ind);
    // frame[89999].print();
    // int ix =0;
    // for(int j=0;j<N;++j){
    //   for(int k=0;k<N;++k){
    //     window[ix]->printall();
    //     dprint("\n");
    //     ++ix;
    //   }
    //   // dprint("\n");
    // }
    // dprint("\n");
    // verbose = 0;

    // dprint("(%d %d %d %d);",x,y,z,ind);
    // if(iteration >= 17670)verbose = 1;

    // if(++iteration%100 == 0){
    //   printf("Done: %d %.8f\n", iteration, (float)((float)iteration/total));
    //   // ++ind;
    //   // continue;
    // }
    // dprint("offset: %d, %d\n",off_x,off_y);

    // dprint("calculating median\n");
    // verbose = 0;
    // printf("Med\n");
    // printf("median of:\n");
    // window.print();
    out[ind] = window.median();
    // printf(" = %d\n", out[ind]);
    // if(iteration > 9600)verbose = 1;

    ind += dir_x;
    x   += dir_x;

    if(x == fmax[0]-R || x == R-1){
      y += dir_y;
      if(y == fmax[1]-R || y == R-1){
        ++z;
        if(z == fmax[2]-R){
          // wrap x, y, z;
          break;
        }else{
          // printf("Y\n");
          // wrap x, y; slide along z and prepare for next line.
          ind = ind + (fmax[0]*fmax[1]) - dir_x;
          x  -= dir_x;
          y  -= dir_y;

          dir_y = (dir_y==1)?-1:1;
          dir_x = (dir_x==1)?-1:1;

          ++depth_z;

          // move window down by updating each
          // of our trees.
          int to_ind, to_xx=off_x, to_yy=off_y, to_zz = off_z;
          int fr_ind=ind-R-(R*dims[0])+(R*dims[0]*dims[1]), fr_xx=0, fr_yy=0;
          
          if(to_xx < 0)to_xx += N;
          if(to_xx >=N)to_xx -= N;
          if(to_yy < 0)to_yy += N;
          if(to_yy >=N)to_yy -= N;
          
          to_ind = to_xx + N*to_yy + N*N*to_zz;

          int ii = 0;
          while(ii < N*N){
            // printf("yo");
            // dprint("(%d, %d) -> %d\n", to_xx, to_yy, to_ind);
            // push the tree down to the next level
            // dprint("add data[%d] to window[%d]\n", fr_ind, to_ind);
            // dprint("set window[%d][%d] = data[%d]\n", to_ind, window[to_ind]->off, fr_ind);
            // vvvv
            // dprint("ff");
            window.replace(to_ind, data[fr_ind]);

            // update indices
            ++to_ind;
            ++to_xx;
            ++fr_ind;
            ++fr_xx;
            if(fr_xx == N){
              // go to next line of data array
              fr_ind += dims[0]-fr_xx;
              fr_xx   = 0;
              ++fr_yy;

              to_ind += N;
              ++to_yy;
              // wrap with offset along y axis.
              if(fr_yy == N){
                fr_yy = 0;
                fr_ind -= N*(N-1);
              }
            }
            // wrap with offset along x and y axes.
            if(to_xx==N){
              to_ind -= to_xx;
              to_xx = 0;
            }
            if(to_yy==N){
              to_ind -= to_yy*(N);
              to_yy = 0;
            }
            // dprint("done\n");
            ++ii;
          }

          // increment z offset.
          ++off_z;
          if(off_z==N)off_z=0;
          if(off_z<0)off_z+=N;
        }
      }else{
        // printf("X %d\n",dir_y);
        // wrap x; slide along y and prepare for next line.
        ind = ind + ((dir_y==1)?fmax[0]:-fmax[0]) - dir_x;
        x  -= dir_x;
        dir_x = (dir_x==1)?-1:1;

        // simple y stutterstep. take the next N
        // trees along the Y axis of frame and put
        // them in window, being aware of the x and
        // y offsets and y_dir.

        int fr_xx=0,     fr_zz=0     ;
        int to_xx=off_x, to_zz=off_z, to_yy= off_y+((dir_y==1)?(0):(-1)); // << todo: replace with (1):(0)?

        if(to_yy<0)to_yy+=N;
        if(to_yy>=N)to_yy-=N;

        int fr_ind = (win_y  + ((dir_y==1)?(N):(-1)) )*fsiz[1] + win_x + depth_z*dims[0]*dims[1];
        int to_ind = to_yy*N + off_x + off_z*N*N;

        // printf("offs -- %d %d %d\n", off_y, to_yy, to_ind);

        // printf("to: %d %d %d %d\n", off_x, off_y, off_z, to_ind);


        if(to_ind>=N*N*N)to_ind -= N*N;
        if(to_ind<0)to_ind += N*N;

        for(int i=0;i<N*N;++i){
          // dprint("set window[%d] = frame[%d]\n",to_ind, fr_ind);
          // printf("%d %d %d\n", to_xx, to_zz, to_ind);
          window.replace(to_ind, data[fr_ind]);

          ++fr_xx;
          ++to_xx;
          fr_ind += 1;
          to_ind += 1;

          // wrap data array input.
          if(fr_xx == N){
            to_ind += N*N;
            fr_ind -= N;
            fr_ind += dims[0]*dims[1];
            fr_xx = 0;
            to_zz += 1;
            ++fr_zz;
          }
          // wrap offsets.
          if(to_xx >= N){
            to_xx = 0;
            to_ind -= N;
          }
          if(to_zz >= N){
            to_zz = 0;
            to_ind -= N*N*N;
          }
        }

        // dprint("xx\n");
        off_y += dir_y;
        if(off_y==N)off_y=0;
        if(off_y<0)off_y+=N;
        win_y += dir_y;
      }
    }else{
      // printf("C %d\n",dir_x);

      int fr_yy=0, fr_zz=0;
      int to_xx=off_x + ((dir_x==1)?(0):(-1)), to_yy=off_y, to_zz=off_z;
      if(to_xx<0)to_xx+=N;
      if(to_xx>=N)to_xx-=N;

      int fr_ind = win_y*fsiz[0] + win_x + ((dir_x==1)?(N):(-1)) + depth_z*dims[0]*dims[1];
      int to_ind = off_y*N       + to_xx + off_z*N*N;

      // printf("ind %d %d %d %d\n", off_x,off_y,off_z,to_ind);

      if(to_ind>=N*N*N)to_ind -= N*N;
      if(to_ind<0)to_ind += N*N;

      for(int i=0;i<N*N;++i){
        // printf("%d %d %d\n", to_yy, to_zz, to_ind);
        // printf("set window[%d] = frame[%d]\n",to_ind, fr_ind);
        window.replace(to_ind, data[fr_ind]);
        ++fr_yy;
        ++to_yy;
        fr_ind += dims[0];
        to_ind += N;

        // wrap data array input.
        if(fr_yy == N){
          to_ind += N*N;
          to_zz += 1;

          fr_ind -= dims[0]*N;
          fr_ind += dims[0]*dims[1];
          fr_yy = 0;
          ++fr_zz;
        }
        // wrap offsets.
        if(to_yy >= N){
          to_yy = 0;
          to_ind -= N*N;
        }
        if(to_zz >= N){
          to_zz = 0;
          to_ind -= N*N*N;
        }
      }

      off_x += dir_x;
      if(off_x==N)off_x=0;
      if(off_x<0)off_x+=N;

      win_x += dir_x;
      // dprint(".\n");
    }
    // printf("P\n");
    // dprint("VV");
  }
  // dprint("\n");
  // free(frame);
}
