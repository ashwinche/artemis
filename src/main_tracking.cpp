#include <SFML/Graphics.hpp>
#include <teem/nrrd.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "experiment.h"
#include "view.h"
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

  sf::Clock clock;
  // sf::Text text;

  bool keys[1024];
  bool clicked[3];
  ivec2 mousepos;

  float scale  = 0;
  Game(){ }

  // "assert" a boolean with a char* message
  void asserts(bool b, const char *message){
    if(!b){
      fprintf(stderr,"ASSERT: %s\n", message);
      ::exit(0);
    }
  }

  // Load font, create and initialize window.
  void initUI(){
    // asserts(font.loadFromFile("../rsc/CallingCode-Regular.ttf"), "loading font");
    // text.setFont(font);
    // text.setString("Loading...");
    // text.setCharacterSize(16);
    // text.setFillColor(sf::Color::White);
    // text.setStyle(sf::Text::Bold);

    window = new sf::RenderWindow(sf::VideoMode(800, 600), "ACellTracking");
    window->setFramerateLimit(30);

    sfview = window->getDefaultView();
  }

  // save the state of the CAMERA
  void save(){
    FILE *file = fopen("save.artemis","wb");
    fwrite(&(view->camera),sizeof(view->camera),1,file);
    fclose(file);
  }

  // load the state of the CAMERA
  void load(){
    FILE *file = fopen("save.artemis","rb");
    if(!file)return;
    int r = fread(&(view->camera),sizeof(view->camera),1,file);
    fclose(file);
  }

  // save the state of the CAMERA and then close the window.
  void exit(){
    save();
    window->close();
  }

  // init camera, create new experiment, init keypresses
  // initialize VIEW (renderer) for a particular EXPERIMENT and computation PIPELINE.

  void init(int argc, char** argv){
    printf("init\n");
    view->camera.set(vec3(-4,3,6), vec3(1,0,-0.33), vec3(0,0,1));

    if(argc == 4){
      experiment = new ArExperiment(argv[1], atoi(argv[2]), atoi(argv[3]), 2);
    }
    else{
      printf("%d", argc);
      printf("USAGE: ./tracking [nrrd_path] [low] [high]\n");
      ::exit(0);
    }

    view->setvolume(experiment->get(0));
    for(int i=0;i<1024;i++)keys[i]=false;
    for(int i=0;i<3;i++)clicked[i]=false;

    // pipeline->process(reprmode.timestep,reprmode.timestep+3);

    load();


    // view->setvolume(pipeline->repr(reprmode));
    // view->setgeometry(pipeline->reprgeometry(reprmode));
    view->touch();
    // printf("initted\n");
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
        view->camera.flat.projmode = '+';
      }
      else if(view->camera.flat.projmode == '+'){
        view->camera.flat.projmode = '_';
      }
      else{
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
    // text.setString(
    //     "(rendered " + to_string(renderframenum) +" frames, " + to_string(ms) + "ms)\n" + 
    //     "gamma:   " + to_string(view->get_gamma()) + "\n" + 
    //     "falloff: " + to_string(view->get_falloff()) + "\n"
    //   );

    window->clear(sf::Color(10,10,10));
    view->render_to(window);


    // window->draw(text);
    window->display();
  }
  int run(int argc, char** argv){
    view = new View(400,400);
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