#include "particle.h"
#include <ctime>
#include <cmath>


Particle::Particle() {
    _radius = 3;
    _mass = 10;
    _surfaceTensionStrength = 20.0;

    _previous_time = 0;

    _rgb = {255, 255, 255};

    _position = {100, 100};
    _velocity = {0.0, 0.0};
    _acceleration = {0.0, 0.0};

    _isColliding = false;
}

Particle::Particle (float radius, float mass, float surfaceTensionStrength, float previous_time, std::vector<int> rgb, std::vector<float> position, std::vector<float> velocity, std::vector<float> acceleration){
    _radius = radius;
    _mass = mass;
    _surfaceTensionStrength = surfaceTensionStrength;

    _previous_time = previous_time;

    _rgb = rgb;

    _position = position;
    _velocity = velocity;
    _acceleration = acceleration;

    _isColliding = false;
}

Particle::Particle(Particle& otherParticle) {
    _radius = otherParticle.getRadius();
    _mass = otherParticle.getMass();

    _previous_time = otherParticle.getPreviousTime();

    _rgb = otherParticle.getColor();

    _position = otherParticle.getPosition();
    _velocity = otherParticle.getVelocity();
    _acceleration = otherParticle.getAcceleration();

    _isColliding = false;
}

std::vector<float> Particle::getPosition () const {
    return this->_position;
}

void Particle::setPosition(std::vector<float> position){
    _position = position;
}

float Particle::getRadius() const {
    return this->_radius;
}

void Particle::setRadius(float radius) {
    if (radius > 0) {
        _radius = radius;
    }
}

float Particle::getMass() const {
    return this->_mass;
}

void Particle::setMass(float mass) {
    if (mass > 0) {
        _mass = mass;
    }
}

float Particle::getSurfaceTensionStrength() const {
  return this->_surfaceTensionStrength;
}

std::vector<int> Particle::getColor() const {
    return this->_rgb;
}

void Particle::setColor(std::vector<int> color) {
    _rgb = color;
}

std::vector<float> Particle::getVelocity() const {
    return this->_velocity;
}

void Particle::setVelocity(std::vector<float> velocity) {
    _velocity = velocity;
}

std::vector<float> Particle::getAcceleration() {
    return this->_acceleration;
}

void Particle::setAcceleration(std::vector<float> acc){
    _acceleration = acc;
}

float Particle::getPreviousTime() {
    return this->_previous_time;
}

void Particle::setPreviousTime(float time) {
    _previous_time = time;
}

bool Particle::getIsColliding() const {
    return this->_isColliding;
}

void Particle::setIsColliding(bool isColliding) {
    _isColliding = isColliding;
}

void Particle::updateParticle(float current_time) {

    float delta_position_x;
    float delta_position_y;

    float elapsed_time = current_time - _previous_time;

    _velocity[0] += _acceleration[0] * elapsed_time;
    _velocity[1] += _acceleration[1] * elapsed_time;

    delta_position_x = _acceleration[0] * elapsed_time + _velocity[0];
    delta_position_y = _acceleration[1] * elapsed_time + _velocity[1];

    _position[0] += delta_position_x;
    _position[1] += delta_position_y;

    _previous_time = current_time;
}
