#include "particle_sim.h"
#include <cmath>
#include <vector>

ParticleSim::ParticleSim(float gravity, float interaction_radius, std::vector<Particle*> particles, int width, int height){
    _gravity = gravity;
    _interaction_radius = interaction_radius;
    _particles = particles;

    _width = width;
    _height = height;
}

void ParticleSim::addParticle(Particle* particle) {
    _particles.push_back(particle);
}

void ParticleSim::updateParticles(float current_time) {
    double pi = M_PI;

    _collidingPairs.clear();

    for (Particle* particle1 : _particles) {

        float Fx_total = 0.0;
        float Fy_total = 0.0;

        for (Particle* particle2 : _particles) {
          Fx_total += _calculateSurfaceTension(particle1, particle2)[0];
          Fy_total += _calculateSurfaceTension(particle1, particle2)[1];
        }

        float friction_force_x = 0.0;
        float friction_force_y = 0.0;
        /*if (_isFrictionEnabled) {*/
        /*    // stokes law for force of drag on a sphere*/
        /*    friction_force_x = 6 * M_PI * particle1->getRadius() * _viscosityOfMedmium * particle1->getVelocity()[0];*/
        /*    friction_force_y = 6 * M_PI * particle1->getRadius() * _viscosityOfMedmium * particle1->getVelocity()[1];*/
        /*}*/

        Fx_total -= friction_force_x;
        Fy_total -= friction_force_y;

        // check for collisions with walls
        // if there is a collision with a wall, negate the force in that direction to account for the normal force
        if (particle1->getPosition()[0] <= 0) {                       // left wall collision
          if (Fx_total < 0.0) {
            Fx_total = 0.0;
          }
        } else if (particle1->getPosition()[0] >= _width) {           // right wall collision
          if (Fx_total > 0.0) {
            Fx_total = 0.0;
          }
        } else if (particle1->getPosition()[1] <= 0) {                // bottom wall collision
          if (Fy_total < 0.0) {
            Fy_total = 0.0;
          }
        } else if (particle1->getPosition()[1] >= _height) {          // top wall collision
          if (Fy_total > 0.0) {
            Fy_total = 0.0;
          }
        }
        
        float a_x = (Fx_total - friction_force_x) / particle1->getMass();
        float a_y = (Fy_total - friction_force_y) / particle1->getMass();

        std::vector<float> acceleration = {a_x, a_y};

        particle1->setAcceleration(acceleration);
    }

    for (std::pair<Particle*, Particle*> collidingPair : _collidingPairs){
      _collision(collidingPair.first, collidingPair.second);
    }

    for (int i = 0; i < _particles.size(); i++){
        _particles[i]->updateParticle(current_time);
    }
}

void ParticleSim::setGravity(const double gravity) {
  _gravity = gravity;
}

const double ParticleSim::getGravity() {
  return _gravity;
}

std::vector<Particle*> ParticleSim::getParticles() const {
    return this->_particles;
}

float ParticleSim::_delta_x(const Particle* particle1, const Particle* particle2) {
    return particle2->getPosition()[0] - particle1->getPosition()[0];
}

float ParticleSim::_delta_y(const Particle* particle1, const Particle* particle2) {
    return particle2->getPosition()[1] - particle1->getPosition()[1];
}

float ParticleSim::_euclideanDistance(const Particle* particle1, const Particle* particle2) {
    return sqrt(pow(_delta_x(particle1, particle2), 2) + pow(_delta_y(particle1, particle2), 2));
}

float ParticleSim::_theta(const Particle* particle1, const Particle* particle2) {
    return atan2(_delta_y(particle1, particle2), _delta_x(particle1, particle2));
}

float ParticleSim::_xComponent(const float value, const float thetaRad) {
    return value * float(cos(thetaRad));
}

float ParticleSim::_yComponent(const float value, const float thetaRad) {
    return value * float(sin(thetaRad));
}

std::vector<float> ParticleSim::_calculateSurfaceTension(Particle* particle1, Particle* particle2) {
  if (particle1 == particle2){    // we do not want to calculate attraction between particles and themselves
    return {0.0, 0.0};
  }

  float euclidean_distance = _euclideanDistance(particle1, particle2);

  float F = 0;

  if (euclidean_distance > _interaction_radius) {   // only calculate attraction for particles that are close enough to interact
    return {0.0, 0.0};
  }

  if (euclidean_distance > (particle1->getRadius() + particle2->getRadius())){        // only calculate force of attraction if particles are not touching
    F = (particle1->getSurfaceTensionStrength() * particle2->getSurfaceTensionStrength()) / pow(euclidean_distance, 2);

    float theta = _theta(particle1, particle2);
    float Fx = _xComponent(F, theta);
    float Fy = _yComponent(F, theta);

    if (particle1->getIsColliding() && particle2->getIsColliding()) {
      particle1->setIsColliding(false);
      particle2->setIsColliding(false);
    }

    return {Fx, Fy};
  } else {    // collision
     if (!particle1->getIsColliding() && !particle2->getIsColliding()) {
      std::pair<Particle*, Particle*> collidingPair = {particle1, particle2};
      _collidingPairs.insert(collidingPair);
    }
    particle1->setIsColliding(true);
    particle2->setIsColliding(true);
    
    return {0.0, 0.0};
  }
}

void ParticleSim::_collision(Particle * particle1, Particle * particle2) {
  float particle1_mass = particle1->getMass();
  float particle2_mass = particle2->getMass();

  float vel1_ix, vel1_iy, vel2_ix, vel2_iy;
  vel1_ix = particle1->getVelocity()[0];
  vel1_iy = particle1->getVelocity()[1];
  vel2_ix = particle2->getVelocity()[0];
  vel2_iy = particle2->getVelocity()[1];

  float vel1_fx, vel1_fy, vel2_fx, vel2_fy;

  vel2_fx = (2 * particle1_mass * vel1_ix + vel2_ix * (particle2_mass - particle1_mass)) / (particle2_mass + particle1_mass);
  vel2_fy = (2 * particle1_mass * vel1_iy + vel2_iy * (particle2_mass - particle1_mass)) / (particle2_mass + particle1_mass);

  vel1_fx = vel2_ix + vel2_fx - vel1_ix;
  vel1_fy = vel2_iy + vel2_fy - vel1_iy;

  std::vector<float> vel1_f, vel2_f; 
  vel1_f = {vel1_fx, vel1_fy};
  vel2_f = {vel2_fx, vel2_fy};

  particle1->setVelocity(vel1_f);
  particle2->setVelocity(vel2_f);
}

std::vector<float> _calculateBuoyantFoce(Particle * particle) {
  
  return {1.0, 1.0};
}

void ParticleSim::reset() {
    _particles.clear();
}

ParticleSim::~ParticleSim(){
    for (Particle* particle : _particles) {
        delete particle;
    }
}
