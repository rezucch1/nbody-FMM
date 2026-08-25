#pragma once

class Morton{
  public:
    static unsigned int encode(const unsigned int* cartesian_coordinates, unsigned int dim);
    static void decode(unsigned int* cartesian_coordinates, unsigned int dim, unsigned code);
};