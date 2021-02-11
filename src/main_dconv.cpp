#include <SFML/Graphics.hpp>
#include <teem/nrrd.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "experiment.h"
#include "view.h"
#include "pipeline.h"
#include "synth.h"
#include <LBFGS.h>


/**
 * Main class which creates and encapsulates the platform
**/
class Game{
public:
  View *view;
  ArExperiment *experiment;
  sf::RenderWindow *window;
  sf::View sfview;
  ArPipeline *pipeline;

  sf::Clock clock;
  sf::Font font;
  sf::Text text;

  ReprMode reprmode;
  bool keys[1024];
  bool clicked[3];
  ivec2 mousepos;

  float scale  = 0;
  Game() : reprmode("plain"){ }

  // "assert" a boolean with a char* message
  void asserts(bool b, const char *message){
    if(!b){
      fprintf(stderr,"ASSERT: %s\n", message);
      ::exit(0);
    }
  }

  // Load font, create and initialize window.
  void initUI(){
    asserts(font.loadFromFile("../rsc/CallingCode-Regular.ttf"), "loading font");
    text.setFont(font);
    text.setString("Loading...");
    text.setCharacterSize(16);
    text.setFillColor(sf::Color::White);
    text.setStyle(sf::Text::Bold);

    window = new sf::RenderWindow(sf::VideoMode(800, 600), "ACellTracking");
    window->setFramerateLimit(30);

    sfview = window->getDefaultView();
  }

  // save the state of the CAMERA
  void save(){
    FILE *file = fopen("../rsc/save.artemis","wb");
    fwrite(&(view->camera),sizeof(view->camera),1,file);
    fclose(file);

    // pipeline->save();
  }

  // load the state of the CAMERA
  void load(){
    FILE *file = fopen("../rsc/save.artemis","rb");
    int r = fread(&(view->camera),sizeof(view->camera),1,file);
    fclose(file);

    // pipeline->load();
  }

  // save the state of the CAMERA and then close the window.
  void exit(){
    save();
    window->close();
  }

  // init camera, create new experiment, init keypresses
  // initialize VIEW (renderer) for a particular EXPERIMENT and computation PIPELINE.

  void init(int argc, char** argv){
    char *volumepath;
    // synth();
    view->camera.set(vec3(-4,3,6), vec3(1,0,-0.33), vec3(0,0,1));

    // experiment = new ArExperiment("/media/ashwin/UBUNTU 18_0/data/???.nrrd",0,9,1);
    // experiment = new ArExperiment("/home/ashwin/data/17-05-01/1??.nrrd",0,9,1);
    // experiment = new ArExperiment("/home/ashwin/data/16-05-05/???.nrrd",0,99,1);

    bool should_track = false;
    if(argc == 4){
      experiment = new ArExperiment(argv[1], atoi(argv[2]), atoi(argv[3]), 2);
      volumepath = argv[1];
      should_track = false;
    }
    else{
      printf("USAGE: ./game [nrrd_path] [min] [max]\n");
    }

    pipeline = new ArPipeline(experiment);

    view->setvolume(pipeline->repr(reprmode));
    for(int i=0;i<1024;i++)keys[i]=false;
    for(int i=0;i<3;i++)clicked[i]=false;

    load(); 

    reprmode.name = "plain";
    reprmode.geom = "paths";

    view->setvolume(pipeline->repr(reprmode));
    view->setgeometry(pipeline->reprgeometry(reprmode));
    view->touch();
  }

  // low-level function to just handle events (at the window level)
  // does not do any processing.
  void handle_events(){
    sf::Event event;
    while (window->pollEvent(event)){
      if (event.type == sf::Event::Closed){
        exit();
      }
      if (event.type == sf::Event::KeyPressed){
        if(event.key.code >= 0){
          keys[event.key.code] = true;
        }
      }
      if (event.type == sf::Event::KeyReleased){
        if(event.key.code >= 0){
          keys[event.key.code] = false;
        }
      }
      if (event.type == sf::Event::Resized){
        window->setView(sfview = sf::View(sf::FloatRect(0,0,window->getSize().x, window->getSize().y)));
      }
      if (event.type == sf::Event::MouseMoved){
        mousepos.x = event.mouseMove.x;
        mousepos.y = event.mouseMove.y;
      }
      if (event.type == sf::Event::MouseButtonPressed){
        if (event.mouseButton.button == sf::Mouse::Left){
          clicked[0] = true;
        }
        if (event.mouseButton.button == sf::Mouse::Right){
          clicked[1] = true;
        }
      }
      if (event.type == sf::Event::MouseButtonReleased){
        if (event.mouseButton.button == sf::Mouse::Left){
          clicked[0] = false;
        }
        if (event.mouseButton.button == sf::Mouse::Right){
          clicked[1] = false;
        }
      }
    }
  }

  // interface between user (keyboard, mouse) and the volume computation pipeline.
  void check_keys(){
    using glm::vec3;
    float speed = 0.1f;

    // exit
    if(keys[sf::Keyboard::F2]){
      exit();
    }
    
    // camera speed
    if(keys[sf::Keyboard::LShift]){
      speed *= 10.f;
    }
    if(keys[sf::Keyboard::LControl]){
      speed *= 0.1f;
    }

    // camera movement
    if(keys[sf::Keyboard::W]){
      view->camera.drawflat = false;
      view->move3D(vec3(0,0,speed));
    }
    if(keys[sf::Keyboard::S]){
      view->camera.drawflat = false;
      view->move3D(vec3(0,0,-speed));
    }
    if(keys[sf::Keyboard::A]){
      view->camera.drawflat = false;
      view->move3D(vec3(-speed,0,0));
    }
    if(keys[sf::Keyboard::D]){
      view->camera.drawflat = false;
      view->move3D(vec3(speed,0,0));
    }
    if(keys[sf::Keyboard::R]){
      view->camera.drawflat = false;
      view->move3D(vec3(0,speed,0));
    }
    if(keys[sf::Keyboard::F]){
      view->camera.drawflat = false;
      view->move3D(vec3(0,-speed,0));
    }
    if(keys[sf::Keyboard::Left]){
      view->camera.drawflat = false;
      view->rotateH(0.1f);
    }
    if(keys[sf::Keyboard::Right]){
      view->camera.drawflat = false;
      view->rotateH(-0.1f);
    }
    if(keys[sf::Keyboard::Down]){
      view->camera.drawflat = false;
      view->rotateV(-0.1f);
    }
    if(keys[sf::Keyboard::Up]){
      view->camera.drawflat = false;
      view->rotateV(0.1f);
    }

    // gamma
    if(keys[sf::Keyboard::Dash]){
      view->step_gamma(1.1f);
      view->touch();
    }
    if(keys[sf::Keyboard::Equal]){
      view->step_gamma(1/1.1f);
      view->touch();
    }

    // light falloff
    if(keys[sf::Keyboard::LBracket]){
      view->step_falloff(1.1f);
      view->touch();
    }
    if(keys[sf::Keyboard::RBracket]){
      view->step_falloff(1/1.1f);
      view->touch();
    }

    // 2D view
    if(keys[sf::Keyboard::M]){
      view->camera.flat.slice += 0.04*speed;
      view->camera.drawflat = true;
      view->touch();
    }
    if(keys[sf::Keyboard::N]){
      view->camera.flat.slice -= 0.04*speed;
      view->camera.drawflat = true;
      view->touch();
    }
    if(keys[sf::Keyboard::Comma]){
      keys[sf::Keyboard::Comma] = false;
      if(view->camera.flat.projmode == 'M'){
        view->camera.flat.projmode = '_';
      }
      else if(view->camera.flat.projmode == '_'){
        view->camera.flat.projmode = '2';
      }else{
        view->camera.flat.projmode = 'M';        
      }
      view->touch();
    }
    if(keys[sf::Keyboard::Period]){
      keys[sf::Keyboard::Period] = false;
      if(view->camera.flat.projaxis == 'x'){
        view->camera.flat.projaxis = 'y';
      }
      else if(view->camera.flat.projaxis == 'y'){
        view->camera.flat.projaxis = 'z';
      }
      else{
        view->camera.flat.projaxis = 'x';
      }
      view->touch();
    }

    // move timestep
    if(keys[sf::Keyboard::O]){
      if(keys[sf::Keyboard::LShift]){
        reprmode.timestep -= 15.f;
      }else if(keys[sf::Keyboard::LControl]){
        reprmode.timestep -= 0.1f;
      }
      else reprmode.timestep -= 1.f;
      view->setvolume(pipeline->repr(reprmode));
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::P]){
      if(keys[sf::Keyboard::LShift]){
        reprmode.timestep += 15.f;
      }else if(keys[sf::Keyboard::LControl]){
        reprmode.timestep += 0.1f;
      }
      else reprmode.timestep += 1.f;
      view->setvolume(pipeline->repr(reprmode));
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }

    // representation mode (for 3D volume rendering)
    if(keys[sf::Keyboard::Num1]){
      reprmode.name = "plain";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Num2]){
      reprmode.name = "dconv";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Num3]){
      reprmode.name = "blobs_all";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Num4]){
      reprmode.name = "diff_blobs";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Num5]){
      reprmode.name = "gaussian";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Num6]){
      reprmode.name = "laplacian";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Num7]){
      reprmode.name = "hessian";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Num8]){
      reprmode.name = "masked";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Num9]){
      reprmode.name = "levelsets";
      long t1 = time(NULL); 
      view->setvolume(pipeline->repr(reprmode));
      printf("elapsed = %ld %ld %ld", (time(NULL)-t1)/3600, (time(NULL)-t1)/60, (time(NULL)-t1));
      view->touch();
    }
    if(keys[sf::Keyboard::Num0]){
      reprmode.name = "show_blobs";
      view->setvolume(pipeline->repr(reprmode));
      view->touch();
    }
    // if(keys[sf::Keyboard::Num7]){
    //   reprmode.name = "blobs_succs";
    //   view->setvolume(pipeline->repr(reprmode));
    //   view->touch();
    // }

    // set scale of blobs to render
    if(keys[sf::Keyboard::I]){
      reprmode = pipeline->repr_coarser(reprmode);
      view->setvolume(pipeline->repr(reprmode));
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::K]){
      reprmode = pipeline->repr_finer(reprmode);
      view->setvolume(pipeline->repr(reprmode));
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }

    // print diagnostic information
    if(keys[sf::Keyboard::B]){
      printf("pos: %.3f %.3f %.3f\n", view->camera.pos.x, view->camera.pos.y, view->camera.pos.z);
      keys[sf::Keyboard::B] = false;
    }

    // render geometry mode
    if(keys[sf::Keyboard::G]){
      if(!strcmp(reprmode.geom, "none")){
        reprmode.geom = "paths";
      }else if(!strcmp(reprmode.geom, "paths")){
        reprmode.geom = "graph";
      }else if(!strcmp(reprmode.geom, "graph")){
        reprmode.geom = "succs";
      }else{
        reprmode.geom = "none";
      }
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
      keys[sf::Keyboard::G] = false;
    }

    if(clicked[0]){
      printf("clicked %d %d\n", mousepos.x, mousepos.y);
      pipeline->repr_highlight(&reprmode, view->get_camera_pos(), view->pixel_to_ray(vec2(mousepos)), keys[sf::Keyboard::LControl], keys[sf::Keyboard::LShift]);
      view->setvolume(pipeline->repr(reprmode, true)); // true = force re-render
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
      clicked[0] = false;
    }
    // if(clicked[1]){
    //   printf("rtclick %d %d\n", mousepos.x, mousepos.y);
    //   pipeline->path_highlight(&reprmode, view->camera.pos*33.f, view->pixel_to_ray(vec2(mousepos)), keys[sf::Keyboard::LControl], keys[sf::Keyboard::LShift]);
    //   // pipeline->repr_highlight(&reprmode, view->camera.pos*33.f, view->pixel_to_ray(vec2(mousepos)), keys[sf::Keyboard::LControl], keys[sf::Keyboard::LShift]);
    //   // view->setvolume(pipeline->repr(reprmode, false)); // true = force re-render
    //   // view->setgeometry(pipeline->reprgeometry(reprmode));
    //   // view->touch();
    //   clicked[1] = false;
    // }

    // process timesteps {n,n+1}
    if(keys[sf::Keyboard::U]){
      pipeline->process(reprmode.timestep,experiment->high);
      // for(int i=0)
      pipeline->saveframe(reprmode.timestep);
      reprmode.name = "plain";
      reprmode.geom = "graph";
      view->setvolume(pipeline->repr(reprmode));
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::J]){
      pipeline->link(reprmode.timestep,reprmode.timestep+30);
      reprmode.name = "blobs";
      reprmode.geom = "succs";
      view->setvolume(pipeline->repr(reprmode));
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::E]){
      // reprmode.name = "blobs";
      pipeline->emit(reprmode, "-blobs.nrrd", reprmode.timestep,reprmode.timestep+100);
    }
    if(keys[sf::Keyboard::T]){
      pipeline->track(reprmode, experiment->filepath);
      view->setvolume(pipeline->repr(reprmode, true));
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::Y]){
      for(int i=experiment->low; i<=experiment->high; i++)
        pipeline->loadframe(i);
      // pipeline->loadframe(reprmode.timestep);

      // pipeline->estimate(reprmode.timestep);
      // pipeline->select_scales(reprmode.timestep);
      // view->setvolume(pipeline->repr(reprmode, true));
      // view->touch();
    }
    if(keys[sf::Keyboard::Tilde]){
      reprmode.highlight.paths = pipeline->paths;
      view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }
    if(keys[sf::Keyboard::L]){
      pipeline->findpaths(1, -1, "quick");
      // else pipeline->findpaths(3, -1, "tgmm");
      // pipeline->process(reprmode.timestep,reprmode.timestep+2);
      // // reprmode.name = "blobs";
      // reprmode.geom = "graph";
      // view->setvolume(pipeline->repr(reprmode));
      // view->setgeometry(pipeline->reprgeometry(reprmode));
      view->touch();
    }
  }
  /* render
       1. debug text;
       2. output from the View;
     and swap buffers
  */
  void renderall(){
    using std::to_string;
    sf::Time elapsed1 = clock.getElapsedTime();
    static int ms = 0;
    static int renderframenum = 0;
    if(view->render()){
      sf::Time elapsed2 = clock.getElapsedTime();
      double time = (elapsed2.asSeconds() - elapsed1.asSeconds());
      ms = time*1000;
      ++renderframenum;
    }
    text.setString(
        "(rendered " + to_string(renderframenum) +" frames, " + to_string(ms) + "ms)\n" + 
        "timestep "+to_string(reprmode.timestep)+"\n"+
        "render mode: " + std::string(reprmode.name) + " - " + std::string(reprmode.geom) + "\n" + 
        "scale:   " + to_string(reprmode.blob.scale) + "\n" + 
        "gamma:   " + to_string(view->get_gamma()) + "\n" + 
        "falloff: " + to_string(view->get_falloff()) + "\n" + 
        ((pipeline->get(reprmode.timestep).complete)?"processed":"")
      );

    window->clear(sf::Color(10,10,10));
    view->render_to(window);


    // window->draw(text);
    window->display();
  }
  int run(int argc, char** argv){
    view = new View(1024,1024);
    init(argc, argv);
    initUI();

    while (window->isOpen()){
      handle_events();
      check_keys();

      renderall();
    }  
  }
};
int main(int argc, char** argv){
  setvbuf(stdout, NULL, _IONBF, 0);  
  Game game;
  return game.run(argc, argv);
}