#include "view.h"

Camera::Camera(){
  lhorz = 1;
  lvert = 1;
  yaw   = 0;
  pitch = 0;

  pos   = vec3();
  look  = vec3();
  up    = vec3();
  right = vec3();
  sky   = vec3();
  
  worldToScreen = mat4();
  drawflat   = false;

  flat.slice = 0;
  flat.projmode   = '_';
  flat.projaxis   = 'z';
}
void Camera::set(vec3 pos, vec3 look, vec3 sky){
  this->pos = pos;
  this->look = normalize(look);
  this->right = normalize(cross(look, sky));
  this->up = cross(right, look);
  this->sky = sky;
}
line3 Camera::to_screen(line3 l, ivec2 screen){
  if(drawflat && flat.projaxis == 'y'){
    line3 o(l.a, l.b);
    o.a.x = l.a.z * 33.f * screen.x/flat.dims.x;
    o.b.x = l.b.z * 33.f * screen.x/flat.dims.x;
    o.a.y = l.a.x * 33.f * screen.y/flat.dims.y;
    o.b.y = l.b.x * 33.f * screen.y/flat.dims.y;
    o.a.z = 1.f;
    o.b.z = 1.f;
    return o;
  }
  if(drawflat && flat.projaxis == 'z'){
    line3 o(l.a, l.b);
    o.a.x = l.a.y * 33.f * screen.x/flat.dims.x;
    o.b.x = l.b.y * 33.f * screen.x/flat.dims.x;
    o.a.y = l.a.x * 33.f * screen.y/flat.dims.y;
    o.b.y = l.b.x * 33.f * screen.y/flat.dims.y;
    o.a.z = 1.f;
    o.b.z = 1.f;
    return o;
  }
  if(drawflat && flat.projaxis == 'x'){
    line3 o(l.a, l.b);
    o.a.x = l.a.y * 33.f * screen.x/flat.dims.x;
    o.b.x = l.b.y * 33.f * screen.x/flat.dims.x;
    o.a.y = l.a.z * 33.f * screen.y/flat.dims.z;
    o.b.y = l.b.z * 33.f * screen.y/flat.dims.z;
    o.a.z = 1.f;
    o.b.z = 1.f;
    return o;
  }
  l.a -= pos;
  l.b -= pos;

  // components in camera-space.
  vec3 ca(dot(l.a,right), dot(l.a,-up), dot(l.a,look));
  vec3 cb(dot(l.b,right), dot(l.b,-up), dot(l.b,look));

  // clip negative z coordinates.
  if(ca.z <= 0 && cb.z <= 0){
    // printf("bad 3\n");
    return line3(vec3(0,0,0),vec3(0,0,0));
  }else{
    if(ca.z <= 0){
      // printf("bad 1 %.2f %.2f %.2f -> %.2f %.2f %.2f\n", ca.x,ca.y,ca.z, cb.x,cb.y,cb.z);
      vec3 v = cb - ca;            // v.z > 0
      ca = ca + v*((-ca.z/v.z)+0.01f);
      // printf("set %.2f %.2f %.2f\n", ca.x,ca.y,ca.z);
    }
    if(cb.z < 0){
      // printf("bad 2\n");
      vec3 v = ca - cb; // v.z > 0
      cb = cb + v*((-cb.z/v.z)+0.01f);
    }
  }

  // project line onto plane distance=1 from camera.
  ca.x /= ca.z;
  ca.y /= ca.z;

  cb.x /= cb.z;
  cb.y /= cb.z;

  l.a.x = (ca[0]/lhorz)*(screen.x/2.f) + (screen.x/2.f);
  l.a.y = (ca[1]/lvert)*(screen.y/2.f) + (screen.y/2.f);
  l.a.z = ca.z;

  l.b.x = (cb[0]/lhorz)*(screen.x/2.f) + (screen.x/2.f);
  l.b.y = (cb[1]/lvert)*(screen.y/2.f) + (screen.y/2.f);
  l.b.z = cb.z;

  return l;
}
// vec3 Camera::to_screen(vec3 x, ivec2 screen){
//   // vec3 p;
  
//   // x = x-pos;


//   // c0 = c0/c2;
//   // c1 = c1/c2; // project onto plane distance 1 from camera.

//   // p.x = (c0/lhorz)*(screen.x/2.f) + (screen.x/2.f);
//   // p.y = (c1/lvert)*(screen.y/2.f) + (screen.y/2.f);
//   // p.z = c2;

//   // return p;
//   return x;
// }