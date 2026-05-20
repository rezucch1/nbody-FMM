#pragma once
#include <memory>
#include <initializer_list>
#include <ostream>


class Tensor{
  public:
    Tensor();
    Tensor(std::initializer_list<double> init);
    Tensor(size_t dim);
    Tensor(const Tensor &other) noexcept;
    Tensor(Tensor &&other) noexcept;
    Tensor &operator=(const Tensor &other);
    Tensor &operator=(const Tensor &&other);

    double &operator[](size_t idx) const;

    Tensor operator+(const Tensor &other) const;
    Tensor operator+(const Tensor &&other) const;

    Tensor operator-() const;
    Tensor operator-(const Tensor &other) const;
    Tensor operator-(const Tensor &&other) const;

    Tensor &operator+=(const Tensor &other);
    Tensor &operator+=(const Tensor &&other);

    Tensor &operator-=(const Tensor &other);
    Tensor &operator-=(const Tensor &&other);

    Tensor &operator*=(double scalar);
    Tensor operator*(double scalar) const;
    Tensor operator/(double scalar) const;
    Tensor &operator/=(double scalar);


    bool operator==(const Tensor &other) const;

    friend std::ostream &operator<<(std::ostream&, const Tensor &);
    

    double squared_norm() const;
    double norm() const;

    size_t dim;
    size_t grade;
    
  protected:

    std::unique_ptr<double[]> data;
};

Tensor operator*(double scalar, const Tensor &vect);
Tensor operator*(double scalar, const Tensor &&vect);

std::ostream &operator<<(std::ostream&, const Tensor&&);
