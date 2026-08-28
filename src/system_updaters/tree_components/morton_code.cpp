/**
 * @file morton_code.cpp
 * @brief Implementation of Morton Z-order curve bit interleaving operations.
 */

#include "morton_code.hpp"
#include <cstring>

size_t Morton::encode(const unsigned int *cartesian_coordinates, unsigned int dim){
  unsigned int code = 0;
  unsigned int temp_cartesian_coordinates[dim];
  std::memcpy(temp_cartesian_coordinates, cartesian_coordinates, sizeof(unsigned int) * dim);
  for (unsigned int i = 0; i < dim; ++i){
    unsigned int digit = i;
    while (temp_cartesian_coordinates[i] > 0){
      code += (temp_cartesian_coordinates[i] % 2) << digit;
      temp_cartesian_coordinates[i] /= 2;
      digit += dim;
    }
  }
  return code;
}

void Morton::decode(unsigned int *cartesian_coordinates, unsigned int dim, size_t code){
  for (unsigned int i = 0; i < dim; ++i)
    cartesian_coordinates[i] = 0;

  unsigned int digit = 0;
  while (code > 0){
    for (unsigned int i = 0; i < dim; ++i){
      cartesian_coordinates[i] += (code % 2) << digit;
      code /= 2;
    }
    ++digit;
  }
}
