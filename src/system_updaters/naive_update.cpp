/**
 * @file naive_update.cpp
 * @brief Implementation of direct O(N^2) gravitational acceleration solver.
 */

#include <iostream>
#include "naive_update.hpp"
#include <cmath>


void NaiveUpdate::update(const std::vector<Particle> &particles) {
  Tensor p_gradient;

  for (const Particle &p : particles) {
    p_gradient = 0 * p.get_position();

    for (const Particle &other : particles) {
      if (&p != &other) {

        Tensor distance = p.get_position() - other.get_position();
        double squared_distance_norm = distance.squared_norm();
        double dim = p.get_position().dim;
        if (dim == 2) {
          p_gradient -= other.get_weight()
                / (squared_distance_norm + 1.0e-4)
                * distance;
}
else if (dim == 3) {
     p_gradient += other.get_weight()
                / ((squared_distance_norm + 1.0e-4) * distance.norm())
                * distance;
}
      }
    }

    p.compute_new_accelaration(p_gradient);
   
  }
}