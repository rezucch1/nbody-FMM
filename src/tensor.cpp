#include "tensor.hpp"
#include <cstring>
#include <iostream>
#include <cmath>

Tensor::Tensor(size_t dim, size_t grade)
  : dim(dim)
  , data(new double[dim])
  {}

Tensor::Tensor(){}

Tensor::Tensor(std::initializer_list<double> init) : Tensor(init.size()){
  std::memcpy(data.get(), init.begin(), sizeof(double) * dim);
}

Tensor::Tensor(const Tensor &other) noexcept : Tensor(other.dim){
  std::memcpy(this->data.get(), other.data.get(), sizeof(double) * dim);
}

Tensor::Tensor(Tensor &&other) noexcept
  : dim(other.dim)
  , data(other.data.release())
  {}

Tensor &Tensor::operator=(const Tensor &other){
  this->dim = other.dim;
  data = std::unique_ptr<double[]>(new double[dim]);
  std::memcpy(data.get(), other.data.get(), sizeof(double) * dim);
  return *this;
}

Tensor &Tensor::operator=(const Tensor &&other){
  return *this = other;
}

double &Tensor::operator[](size_t idx) const{
  if (idx >= dim){
    std::cerr << "You have access outside the dimention of the tensor" << std::endl;
    std::exit(1);
  }

  return data[idx];
}

Tensor &Tensor::operator+=(const Tensor &other){
  if (other.dim != this->dim){
    std::cerr << "Can be summed only tensors with the same dim" << std::endl;
    std::exit(1);
  }

  for (size_t i = 0; i < dim; ++i)
    this->data[i] += other.data[i];

  return *this;
}

Tensor &Tensor::operator+=(const Tensor &&other){
  return *this += other;
}

Tensor Tensor::operator+(const Tensor &other) const{
  Tensor res(*this);
  return res+=other;
}

Tensor Tensor::operator+(const Tensor &&other) const{
  return *this + other;
}

Tensor &Tensor::operator-=(const Tensor &other){
  if (other.dim != this->dim){
    std::cerr << "Can be summed only tensors with the same dim" << std::endl;
    std::exit(1);
  }

  for (size_t i = 0; i < dim; ++i)
    this->data[i] -= other.data[i];

  return *this;
}

Tensor &Tensor::operator-=(const Tensor &&other){
  return *this -= other;
}

Tensor Tensor::operator-(const Tensor &other) const{
  Tensor res(*this);
  return res-=other;
}

Tensor Tensor::operator-(const Tensor &&other) const{
  return *this - other;
}

Tensor Tensor::operator-() const{
  Tensor other(this->dim);
  for (size_t i = 0; i < dim; ++i)
    other.data[i] = -this->data[i];
  return other;
}

Tensor &Tensor::operator*=(double scalar){
  for (size_t i = 0; i < dim; ++i)
    this->data[i] *= scalar;
  return *this;
}

Tensor Tensor::operator*(double scalar) const{
  Tensor other(*this);
  return other *= scalar;
}

Tensor operator*(double scalar, const Tensor &vect){
  return vect*scalar;
}

Tensor operator*(double scalar, const Tensor &&vect){
  return vect*scalar;
}

Tensor &Tensor::operator/=(double scalar)
{
  for (size_t i = 0; i < dim; ++i)
    data[i] /= scalar;
  return *this;
}

Tensor Tensor::operator/(double scalar) const{
  Tensor other(*this);
  return other /= scalar;
}

bool Tensor::operator==(const Tensor &other) const{
  if (this->dim != other.dim)
    return false;
  return memcmp(this->data.get(), other.data.get(), this->dim * sizeof(double)) == 0;
}

std::ostream &operator<<(std::ostream &stream, const Tensor &vect){
  stream << "(";
  if (vect.dim >= 1){
    stream << vect.data[0];
    for (size_t i = 1; i < vect.dim; ++i)
    stream << ", " << vect.data[i];
  }
  return stream << ")";
}

std::ostream &operator<<(std::ostream &stream, const Tensor &&vect){
  return stream << vect;
}

double Tensor::squared_norm() const{
  double res;
  for (size_t i = 0; i < dim; ++i){
    double value = data[i];
    res += value * value;
  }
  return res;
}

double Tensor::norm() const{
  return std::sqrt(squared_norm());
}