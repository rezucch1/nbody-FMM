/**
 * @file particle.cpp
 * @brief Implementation of Particle methods.
 */

#include "particle.hpp"

Particle::Particle(double weight, Tensor &&pos, Tensor &&vel)
	: weight(weight)
	, acceleration(0 * pos)
	, position(pos)
	, velocity(vel)
	{}

const Tensor &Particle::compute_new_accelaration(const Tensor &potential_gradient) const{
  acceleration = - GRAV_CONST * potential_gradient;
  return acceleration;
}

const Tensor &Particle::get_acceleration() const {
  return acceleration;
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