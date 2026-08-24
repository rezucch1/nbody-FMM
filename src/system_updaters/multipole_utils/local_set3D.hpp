#pragma once
#include "local_set.hpp"
template<>

class LocalSet<3> : public LocalSetI{
  public:
    LocalSet(unsigned int L);
    virtual std::complex<double> operator()(unsigned l, int m=0) const override;
    virtual LocalSet<3> &operator+=(const LocalSetI *other) override;
    virtual LocalSet<3> &operator+=(const LocalSet<3> other);
    virtual LocalSet<3> *distribute_parent_with_distance(const Tensor &d) const override;
    void set_elements(unsigned int n, int m, std::complex<double>);
};


extern template class LocalSet<3>;