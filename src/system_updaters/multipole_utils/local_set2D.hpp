#pragma once

#include"local_set.hpp"

template<>
class LocalSet<2> : public LocalSetI{
  public:
    LocalSet(unsigned int L);
    virtual std::complex<double> operator()(unsigned l, int m=0) const override;
    std::complex<double> &get_element(unsigned int i);
    const std::complex<double> &get_element(unsigned int i) const;
    virtual LocalSet<2> &operator+=(const LocalSetI *other) override;
    virtual LocalSet<2> &operator+=(const LocalSet<2> other);
    virtual LocalSet<2> *distribute_parent_with_distance(const Tensor &d) const override;
    virtual Tensor get_gradient(const Tensor &d) const override;
  protected:
    std::vector<std::complex<double>> elements;
};


extern template class LocalSet<2>;
