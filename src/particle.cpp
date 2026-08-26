#include "particle.hpp"

Particle::Particle(double weight, Tensor &&pos, Tensor &&vel)
	: weight(weight)
	, position(pos)
	, velocity(vel)
	{}

const Tensor Particle::get_acceleration(const Tensor &potential_gradiend) const{
  return - GRAV_CONST * potential_gradiend;
}

double Particle::get_weight() const
{
	return weight;
}

const Tensor &Particle::get_position() const{
  return position;
}

const Tensor &Particle::get_velocity() const{
  return velocity;
}