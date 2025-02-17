#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include <cmath>
#include <sstream>
#include <random>
#include <iostream>

#include "particle.h"
#include "particle_sim.h"

using namespace sf;
using namespace std;

const int WIDTH = 1000;
const int HEIGHT = 1000;

const int MENU_WIDTH = 400;
const int MENU_HEIGHT = 800;


mt19937 mt(time(nullptr));

struct testParticle{
  float radius = 4.0;
  float mass = 5.0;
  vector<int> rgb = {255, 0, 0};
  vector<float> position = {int(WIDTH / 2), int(HEIGHT / 2)};
  vector<float> velocity = {0.0, 0.0};
  vector<float> acceleration = {0.0, 0.0};
  float surface_tension_strength = 10;
};

// convert coordinate from between window and sim coordinates (invert y axis)
vector<float> convertCoords(vector<float> coords){
  return {coords[0], HEIGHT - coords[1]};
}

Vector2f convertCoords(Vector2f coords){
  return {coords.x, HEIGHT - coords.y};
}


int main(){

  float gravity = 9.81;
  float interactionRadius = 10;

  vector<Particle*> particles;

  // mouse position in config menu coordinates
  Vector2i mousePositionMenu; 
  // mouse position in particle menu coordinates
  vector<float> mousePosition = {-10.0, -10.0};
  bool spawning_particle = false;

  testParticle particle_struct;

  Font font;
  FileInputStream fontIn;
  fontIn.open("slkscr.ttf");
  font.loadFromStream(fontIn);

  // RENDER WINDOWS 
  // =======================================================================
  RenderWindow particleWindow(VideoMode(WIDTH, HEIGHT), "Particle Sim");
  particleWindow.setFramerateLimit(60);
  // =======================================================================

  bool sim_running = true;
  bool space_pressed = false;

  Text pausedIndicator; 
  pausedIndicator.setFont(font);
  pausedIndicator.setFillColor(Color::Red);
  pausedIndicator.setPosition(Vector2f(5, 5));
  pausedIndicator.setCharacterSize(30);
  pausedIndicator.setString("Paused.");

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> particle_offset_uniform_dist(0, 20);

  Clock renderTime;

  vector<CircleShape> particle_shapes;

  ParticleSim particleSim(gravity, interactionRadius, particles, WIDTH, HEIGHT);

  // spawn a particle with a fixed position at first, until mouse released then unfix position 
  Particle particle = Particle(particle_struct.radius, particle_struct.mass, particle_struct.surface_tension_strength, 0, particle_struct.rgb, convertCoords(mousePosition), particle_struct.velocity, particle_struct.acceleration);

  Particle * particle_pointer = nullptr;
  
  CircleShape particle_shape;
  stringstream info_stream;

  int renderIteration = 0;
  while(particleWindow.isOpen()){

    Time time = renderTime.getElapsedTime();
    float time_seconds = time.asSeconds();

    mousePosition = {float(Mouse::getPosition(particleWindow).x), float(Mouse::getPosition(particleWindow).y)};
    
    // if mouse pressed within bounds of render window and render window has OS focus
    if (Mouse::isButtonPressed(Mouse::Left) && mousePosition[0] < WIDTH && mousePosition[0] >= 0 && mousePosition[1] < HEIGHT && mousePosition[1] >= 0 && particleWindow.hasFocus()){
      if (renderIteration % 5  == 0) {
        particle.setColor(particle_struct.rgb);
        
        float x_offset = float(particle_offset_uniform_dist(gen)) - 10.0;
        float y_offset = float(particle_offset_uniform_dist(gen)) - 10.0;

        vector<float> particle_spawn_pos = convertCoords(mousePosition);
        particle_spawn_pos[0] += x_offset;
        particle_spawn_pos[1] += y_offset;

        particle.setPosition(particle_spawn_pos);
        particle.setPreviousTime(time_seconds);
        particle.setVelocity({0.0, 0.0, 0.0});
        particle_pointer = new Particle(particle);
        particleSim.addParticle(particle_pointer);
      }
    } else {
      particle_pointer = nullptr;
    }

    // CLOSE WINDOWS IF X PRESSED
    // ==========================================================
    Event particleWindowEvent;
    Event menuWindowEvent;

    while(particleWindow.pollEvent(particleWindowEvent)){
      if(particleWindowEvent.type == Event::Closed){
        particleWindow.close();
      }
    }
    // ==========================================================

    // UPDATE PARTICLES AND DISPLAY
    // ==========================================================

    // KEYBOARD EVENTS =========================================
    if (Keyboard::isKeyPressed(Keyboard::Space)){   // space to pause/unpause
      if (!space_pressed) {
        sim_running = !sim_running;
        space_pressed = true;
      }
    } else {
      space_pressed = false;
    }

    if (Keyboard::isKeyPressed(Keyboard::R)) {
      particleSim.reset();
    }

    particleWindow.clear();

    if (sim_running) {
      particleSim.updateParticles(time_seconds);
    } else {
      particleWindow.draw(pausedIndicator);
    }

    particles = particleSim.getParticles();
    for (int i = 0; i < particles.size(); i ++){
      Particle* particle = particles[i];

      Vector2f particle_pos_sim = Vector2f(particle->getPosition()[0], particle->getPosition()[1]);
      Vector2f particle_pos_window = convertCoords(particle_pos_sim);

      particle_pos_window = Vector2f(particle_pos_window.x - particle->getRadius(), particle_pos_window.y - particle->getRadius());

      particle_shape.setRadius(particle->getRadius());
      particle_shape.setFillColor(Color(particle->getColor()[0], particle->getColor()[1], particle->getColor()[2]));
      particle_shape.setPosition(particle_pos_window);

      particleWindow.draw(particle_shape);
    }

    particleWindow.display();
    // ==========================================================

    renderIteration ++;
  }

  return 0;
}
