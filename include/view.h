#ifndef VIEW_H
#define VIEW_H

#include <glm/glm.hpp>
#include <SFML/Graphics.hpp>
#include "colormap.h"
#include "experiment.h"
#include "shapes.h"

using namespace glm;

struct line3{
  vec3 a;
  vec3 b;

  line3(vec3 a, vec3 b): a(a), b(b){}
  vec3 &operator[](int i){
    if(i==0)return a;
    if(i==1)return b;
    fprintf(stderr, "invalid access: line3[%d]\n", i);
    exit(0);
  }
};

class Camera{
public:
  Camera();
  float yaw,pitch;
  float lhorz, lvert;
  vec3 pos;
  vec3 look, up, right;
  vec3 sky;
  mat4 worldToScreen;

  bool drawflat;

  struct{
    float slice;
    vec3 dims;
    char projmode;
    char projaxis;
  }flat;
  struct{
    bool drawplane;
  }volume;

  void setYawPitch(float y, float p);
  void set(vec3 pos, vec3 look, vec3 up);
  line3 to_screen(line3 x, ivec2 screen);
};

/* Describes a View into a particular ArPipeline.
 * We can view a volume either as
    -  a volume-rendered mesh (with a floating camera); or 
    -  a series of 2D slices that can be traversed along the 3rd dimension.
   We can also render planar geometry, as arbitrary line segments in the ArGeometry class.
 */
class View{
public:
  View(int w, int h);
private:
  ArGeometry3D geom;

  sf::VertexArray lines;
  // sf::VertexArray lines  _window;

  sf::Uint8 *texdata;
  int const w, h;
  sf::Sprite  sprite;
  sf::Texture texture;
  int beat;
  int unstable;

  float gamma;
  float falloff;
  float render_scale;

  struct{
    Nrrd *n;
    NrrdAxisInfo *a;
    int w0,w1,w2,w3;
    int    a1,a2,a3;
    float *data;
  }vcache;

  struct{
    int width;
    int height;
  }win;

  void drawflat();
  void draw_geometry();
  void raytrace();
  float qsample(int c, float x, float y, float z);

public:
  Camera camera;
  Colormap colormap;
  float position;

  // colormap
  void step_gamma(float factor);

  // render
  void step_falloff(float factor);
  void touch();
  int  render();

  // multithreading
  static void *t_raytrace(void*);

  // input
  void setvolume(Nrrd *nrrd);           // set the volume to be rendered
  void setgeometry(ArGeometry3D*);      // set the planar geometry to be rendered
  vec3 pixel_to_ray(vec2);
  vec3 get_camera_pos();

  // movement
  void move3D(vec3 v);
  void rotateH(float r);
  void rotateV(float r);

  // output
  void render_to(sf::RenderWindow *window);
  float get_gamma();
  float get_falloff();
};

#endif