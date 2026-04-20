#include "particle.hpp"

Particle::Particle(double mass, Tensor &&pos, Tensor &&vel)
	: mass(mass)
	, position(pos)
	, velocity(vel)
	{}

const double &Particle::get_mass() const{
  return mass;
}

const Tensor &Particle::get_position() const{
  return position;
}

const Tensor &Particle::get_velocity() const{
  return velocity;
}