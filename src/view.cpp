#include <SFML/Graphics.hpp>
#include <teem/nrrd.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <GLFW/glfw3.h>

#include <pthread.h>
#include <stdio.h>


#include "view.h"
#include "colormap.h"

static inline int ftoi(float i){
  return int(i*255.999999);
}
static inline float itof(int i){
  return float(i)/255.999999;
}
static inline float sq(float in){
  return in*in;
}
static inline void put_color(float x, Colormap &colormap, sf::Uint8 *color){
  if(x<0)x=0;
  if(x>=6.f){
    x -= 6.f;
  }else if(x>=4.f){
    x -= 4.f;
    color[0] = ftoi(x);
    color[2] = ftoi(sqrt(x));
    color[1] = ftoi(x*x);
    color[3] = 255;
  }else if(x>=2.f){
    x -= 2.f;
    color[1] = ftoi(x);
    color[0] = ftoi(sqrt(x));
    color[2] = ftoi(x*x);
    color[3] = 255;
  }else{
    color[0] = ftoi(x);
    color[1] = ftoi(sqrt(x));
    color[2] = ftoi(x*x);
    color[3] = 255;
  }
  // if(x>=1)x=1;
  // vec4 col = colormap.colorof(x);
  // color[0] = ftoi(col.x);
  // color[1] = ftoi(col.y);
  // color[2] = ftoi(col.z);
  // color[3] = ftoi(col.w);
  // color[0] = ftoi(x);
  // color[1] = ftoi(sqrt(x));
  // color[2] = ftoi(x*x);
  // color[3] = 255;

}

View::View(int w, int h) : w(w), h(h), colormap(4.5f), gamma(4.5f), falloff(0.001f){
  texture.create(w,h);
  texdata  = new sf::Uint8[w*h*4];
  camera   = Camera();
  unstable = 1;
  beat     = 0;
  camera.volume.drawplane = false;
  
}

float View::get_gamma(){
  return gamma;
}

void View::step_gamma(float factor){
  if(factor<=0)return;
  colormap.destroy();
  gamma *= factor;
  colormap = Colormap(gamma);
}

float View::get_falloff(){
  return falloff;
}

void View::step_falloff(float factor){
  if(factor<=0)return;
  falloff *= factor;
}

void View::setvolume(Nrrd *n){
  // printf("View::setvolume = %p\n", n);
  if(vcache.n == n){
    return;
  }
  NrrdAxisInfo *a = n->axis;

  vcache.n = n;
  vcache.a = n->axis;
  vcache.data = (float*)n->data;
  vcache.w0 = 1;                      // crawl x
  vcache.w1 = a[0].size * vcache.w0;  // crawl y
  vcache.w2 = a[1].size * vcache.w1;  // crawl z
  vcache.w3 = a[2].size * vcache.w2;  // total size of data

  vcache.a1 = a[0].size;
  vcache.a2 = a[1].size;
  vcache.a3 = a[2].size;

  // printf("vcache %d %d %d\n", n->axis[0].size, n->axis[1].size, n->axis[2].size);
  touch();
}
void View::setgeometry(ArGeometry3D* geometry){
  geom = *geometry;
  lines = sf::VertexArray(sf::Lines, geom.lines.size());
}
void View::draw_geometry(){
  float square = min(win.width, win.height);
  float px = (win.width-square)/2.f;
  float py = (win.height-square)/2.f;

  for(int i=0;i<geom.lines.size();i+=2){
    lines[i+0].color = geom.lines_c[i+0];
    lines[i+1].color = geom.lines_c[i+1];

    line3 clipped = camera.to_screen(line3(geom.lines[i]/33.f, geom.lines[i+1]/33.f), ivec2(square, square));


    // vec3 ps0 = camera.to_screen(geom.lines[i+0], ivec2(square, square));
    // vec3 ps1 = camera.to_screen(geom.lines[i+1], ivec2(square, square));

    if(clipped[0].z <= 0 && clipped[1].z <= 0){
      // if both are behind the camera, then discard both.
      lines[i  ].color = sf::Color::Transparent;
      lines[i+1].color = sf::Color::Transparent;
    }

    bool too_far = clipped[1].z > 15 || clipped[0].z > 15;
    too_far = false;
    if(too_far){
      lines[i  ].color = sf::Color::Transparent;
      lines[i+1].color = sf::Color::Transparent;
      lines[i+0].position.x = 0;
      lines[i+0].position.y = 0;
      lines[i+1].position.x = 0;
      lines[i+1].position.y = 0;
      continue;
    }
    lines[i  ].color.a = 100;
    lines[i+1].color.a = 100;
    // else{
    //   if(ps0.z < 0){
    //     vec3 v = ps1 - ps0; // v.z > 0
    //     ps0 = ps0 + v*(-ps0.z/v.z);
    //   }
    //   if(ps1.z < 0){
    //     vec3 v = ps0 - ps1; // v.z > 0
    //     ps1 = ps1 + v*(-ps1.z/v.z);
    //   }
    // }

    // transform to screen space:
    lines[i+0].position.x = clipped[0].x + px;
    lines[i+0].position.y = clipped[0].y + py;
    lines[i+1].position.x = clipped[1].x + px;
    lines[i+1].position.y = clipped[1].y + py;

    // printf("line %.3f %.3f -> %.3f %.3f\n", lines[i].position.x, lines[i].position.y, lines[i+1].position.x, lines[i+1].position.y);
  }
  
  // lines[0].position = camera.to_screen(vec3());sf::Vector2f(0,0);
  // lines[1].position = sf::Vector2f(0.5f,0.5f);

}
float View::qsample(int c, float x, float y, float z){
  int i0 = c;
  int i1 = x;
  int i2 = y;
  int i3 = z;
  if(i1<0 || i2<0 || i3<0 || i1 >= vcache.a1 || i2 >= vcache.a2 || i3 >= vcache.a3){
    return -1.f;
  }
  if(camera.volume.drawplane && !camera.drawflat){
    if(camera.flat.projaxis == 'x'){
      if(x > camera.flat.slice*vcache.a1)return 0;
    }
    if(camera.flat.projaxis == 'y'){
      if(y > camera.flat.slice*vcache.a2)return 0;
    }
    if(camera.flat.projaxis == 'z'){
      if(z > camera.flat.slice*vcache.a3)return 0;
    }
  }
  // printf("sample %d\n", i3*vcache.w3 + i2*vcache.w2 + i1*vcache.w1 + i0);
  return vcache.data[
    i3*vcache.w2 + i2*vcache.w1 + i1*vcache.w0
  ];
}

void View::drawflat(){
  camera.flat.dims = vec3(vcache.a1, vcache.a2, vcache.a3);
  if(camera.flat.slice>=1)camera.flat.slice=1;
  if(camera.flat.slice<0)camera.flat.slice=0;
  float zz=camera.flat.slice;
  int i=0;
  char axis = camera.flat.projaxis;
  // printf("axis=%c\n", camera.flat.projaxis);
  float ww=w, hh=h;
  if(axis == 'x'){
    ww = vcache.a2;
    hh = vcache.a3;
  }
  if(axis == 'y'){
    ww = vcache.a1;
    hh = vcache.a3;
  }
  if(axis == 'z'){
    ww = vcache.a1;
    hh = vcache.a2;
  }
  if(ww != w){
    float scale = w/ww;
    ww *= scale;
    hh *= scale;
  }
  if(hh > h){
    float scale = h/hh;
    ww *= scale;
    hh *= scale;
  }
  float xx=0,yy=0, vv=0;
  for(int x=0;x<w;++x){
    for(int y=0;y<h;++y){
      vv = 0;
      if(x<ww && y<hh){
        if(axis == 'z'){
          xx = float(x)/ww * (float(vcache.a1)-0.00005f);
          yy = float(y)/hh * (float(vcache.a2)-0.00005f);
          zz = camera.flat.slice * (float(vcache.a3)-0.00005f);
        }if(axis == 'x'){
          xx = camera.flat.slice * (float(vcache.a1)-0.00005f);
          yy = float(x)/ww * (float(vcache.a2)-0.00005f);
          zz = float(y)/hh * (float(vcache.a3)-0.00005f);
        }if(axis == 'y'){
          xx = float(x)/ww * (float(vcache.a1)-0.00005f);
          yy = camera.flat.slice * (float(vcache.a2)-0.00005f);
          zz = float(y)/hh * (float(vcache.a3)-0.00005f);
        }
        if(camera.flat.projmode == 'M'){
          float sum = 0;
          for(int i=0;i<vcache.a3;i++){
            if(axis == 'x')sum = max(sum, qsample(0,i,yy,zz));
            if(axis == 'y')sum = max(sum, qsample(0,xx,i,zz));
            if(axis == 'z')sum = max(sum, qsample(0,xx,yy,i));
          }
          vv = sum;
        }else if(camera.flat.projmode == '+'){
          float sum = 0;
          for(int i=0;i<vcache.a3;i++){
            if(axis == 'x')sum += qsample(0,i,yy,zz);
            if(axis == 'y')sum += qsample(0,xx,i,zz);
            if(axis == 'z')sum += qsample(0,xx,yy,i);
          }
          vv = sum / vcache.a3;
        }else if(camera.flat.projmode == '2'){
          // printf("2");
          float sum = 0;
          int n = 0;
          if(axis == 'x'){
            for(int i=max(0,int(xx-30)); i<min(int(xx+30),vcache.a1); i++){
              sum += qsample(0,i,yy,zz);
              ++n;
            }
          }
          if(axis == 'y'){
            for(int i=max(0,int(yy-30)); i<min(int(yy+30),vcache.a2); i++){
              sum += qsample(0,xx,i,zz);
              ++n;
            }
          }
          if(axis == 'z'){
            for(int i=max(0,int(zz-30)); i<min(int(zz+30),vcache.a3); i++){
              sum += qsample(0,xx,yy,i);
              ++n;
            }
          }
          vv = sum / n;
        }else{
          if(axis == 'x')vv = qsample(0,xx,yy,zz);
          if(axis == 'y')vv = qsample(0,xx,yy,zz);
          if(axis == 'z')vv = qsample(0,xx,yy,zz);
        }
      }
     // vv = sum*5.f;
      // vv = qsample(0,zz,xx,yy);
      put_color(vv, colormap, texdata+i);
      i+=4;
    }
  }
  texture.update(texdata);
}
vec3 View::pixel_to_ray(vec2 v){
  
  float square = min(win.width, win.height);
  float px = (win.width-square)/2.f;
  float py = (win.height-square)/2.f;

  vec2 screen((v.x-px)/square, (v.y-py)/square);

  screen.x = (screen.x-0.5f)*2.f;
  screen.y = (screen.y-0.5f)*2.f;

  // vec3 left = -camera.right*camera.lhorz;
  // vec3 up   = camera.up*camera.lvert;

  // vec3 topleft = forward - (right*camera.lhorz) + (up*camera.lvert);

  printf("clicked %.2f %.2f\n",screen.x, screen.y);
  return vec3(camera.look + camera.right*screen.x*camera.lhorz - camera.up*screen.y*camera.lvert);
  // return camera.look;
}
vec3 View::get_camera_pos(){
  return camera.pos*33.f;
}
struct render_frame_info{
  View *view;
  vec3 startp;
  vec3 topleft;
  vec3 dx;
  vec3 dy;
};
struct thread_raytrace_info{
  render_frame_info *fi;
  int py_min;   // render window min y value.
  int py_max;   // render window max y value.
};
void* View::t_raytrace(void* vinfo){
  thread_raytrace_info *info = (thread_raytrace_info*)vinfo;

  vec3 ray;
  vec3 &topleft = info->fi->topleft;
  vec3 &dx      = info->fi->dx;
  vec3 &dy      = info->fi->dy;
  vec3 &startp  = info->fi->startp;

  vec3 p;
  float color_x,color_y,color_z,color_w;
  vec4 probe;

  int w = info->fi->view->w;
  int h = info->fi->view->h;
  int i = 4*w*info->py_min;
  for(int py=info->py_min;py<info->py_max;++py){
    for(int px=0;px<w;++px){
      ray = topleft + float(px)*dx + float(py)*dy;

      ray = ray * 0.033f * 33.f * 0.5f;
      p   = startp;

      color_x = 0;
      color_y = 0;
      color_z = 0;
      color_w = 1;

      if(p.x < 0)p = p + ray*(-p.x/ray.x);
      if(p.y < 0)p = p + ray*(-p.y/ray.y);
      if(p.z < 0)p = p + ray*(-p.z/ray.z);


      if(p.x > info->fi->view->vcache.a1)p = p - ray*(p.x - info->fi->view->vcache.a1)/ray.x;
      if(p.y > info->fi->view->vcache.a2)p = p - ray*(p.y - info->fi->view->vcache.a2)/ray.y;
      if(p.z > info->fi->view->vcache.a3)p = p - ray*(p.z - info->fi->view->vcache.a3)/ray.z;

      while(color_w > 0.01f){
        p += ray;
  
        float v = info->fi->view->qsample(0, p.x, p.y, p.z);
        if(v<0)break;
        // if(v>0)printf("v=%.2f\n",v);

        if(v > 0.001){
          probe = info->fi->view->colormap.colorof(v);
          if(probe.w > 0.01f){
            color_x += probe.x*probe.w*color_w;
            color_y += probe.y*probe.w*color_w;
            color_z += probe.z*probe.w*color_w;
            color_w += -color_w*probe.w;
          }
        }

        color_w -= info->fi->view->falloff; // some light is absorbed or refracted away.
        
        // color_w *= 0.995f;           
      }
      info->fi->view->texdata[i+0] = int(color_x*255.999999);
      info->fi->view->texdata[i+1] = int(color_y*255.999999);
      info->fi->view->texdata[i+2] = int(color_z*255.999999);
      info->fi->view->texdata[i+3] = 255;

      // info->view->texdata[i+0] = 255;
      // info->view->texdata[i+1] = 128;
      // info->view->texdata[i+2] = 0;
      // info->view->texdata[i+3] = 255;
      i+=4;
    }
  }
  return 0;
}
void View::raytrace(){
  vec3 forward = camera.look;
  vec3 right   = camera.right;
  vec3 up      = camera.up;

  vec3 ray = forward;
  vec3 dx  = (right * camera.lhorz)/(float(w)/2.f);
  vec3 dy  = -(up   * camera.lvert)/(float(h)/2.f);
  vec3 topleft = forward - (right*camera.lhorz) + (up*camera.lvert);

  // printf("raytrace, %.2f %.2f %.2f + %.2f %.2f %.2f\n",camera.pos.x, camera.pos.y, camera.pos.z, forward.x,forward.y,forward.z);

  int i=0;

  vec3 p;
  float color_x,color_y,color_z,color_w;
  vec4 probe;
  vec3 startp = camera.pos * 33.f;

  // multithreaded code.
  const int nthreads = 24;
  pthread_t threads[nthreads];
  thread_raytrace_info info[nthreads];
  render_frame_info fi;
  fi.view = this;
  fi.startp = startp;
  fi.topleft = topleft;
  fi.dx = dx;
  fi.dy = dy;
  for(int i=0;i<nthreads;++i){
//     struct render_frame_info{
//   View *view;
//   vec3 startp;
//   vec3 topleft;
//   vec3 dx;
//   vec3 dy;
// };
    info[i].fi = &fi;
    info[i].py_min = (h*i)/nthreads;
    info[i].py_max = (h*(i+1))/nthreads;
    if(pthread_create(threads+i, NULL, t_raytrace, info+i)){
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

  texture.update(texdata);
}
void View::touch(){
  unstable = 2;
}
int View::render(){
  if(unstable<=0)return (unstable = 0);
  if(camera.drawflat){
    drawflat();
  }else{
    raytrace();
  }
  draw_geometry();
  --unstable;
  return 1;
}

void View::move3D(vec3 v){
  touch();
  // printf("move: %.2f %.2f %.2f\n",v.x,v.y,v.z);
  // printf("camera: %.2f %.2f %.2f\n",camera.right.x,camera.right.y,camera.right.z);
  camera.pos += v.x*camera.right + v.y*camera.up + v.z*camera.look;
}
void View::rotateH(float r){
  touch();
  camera.set(camera.pos, camera.look - camera.right*r, camera.sky);
}
void View::rotateV(float r){
  touch();
  float dot = glm::dot(camera.look,camera.sky);
  if((r < 0 && dot>-0.995) || (r > 0 && dot < 0.995)){
    camera.set(camera.pos, camera.look + camera.up*r, camera.sky);
  }
}
void View::render_to(sf::RenderWindow *window){
  static sf::VertexArray quad(sf::TriangleFan, 4);

  sf::Vector2u win_s = window->getSize();

  if(win_s.x != win.width || win_s.y != win.height){
    win.width  = win_s.x;
    win.height = win_s.y;
    draw_geometry();
  }

  float square = min(win_s.x, win_s.y);
  float px = (win_s.x-square)/2.f;
  float py = (win_s.y-square)/2.f;

  quad[0].position = sf::Vector2f(px,py);
  quad[1].position = sf::Vector2f(px+square,py);
  quad[2].position = sf::Vector2f(px+square, py+square);
  quad[3].position = sf::Vector2f(px,py+square);

  quad[0].texCoords = sf::Vector2f(0,0);
  quad[1].texCoords = sf::Vector2f(w,0);
  quad[2].texCoords = sf::Vector2f(w,h);
  quad[3].texCoords = sf::Vector2f(0,h);

  window->draw(quad, &texture);
  // glLineWidth(4.f);
  window->draw(lines);
}