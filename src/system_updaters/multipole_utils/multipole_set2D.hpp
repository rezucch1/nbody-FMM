#pragma once

#include "multipole_set.hpp"

template<>
class MultipoleSet<2> : public MultipoleSetI {
  
    public: 
    explicit MultipoleSet(unsigned int);
    virtual std::complex<double> operator()(unsigned l, int m=0) const override;
    virtual std::unique_ptr<MultipoleSetI> weigh_children_with_distance(const Tensor &d) const override;
    virtual std::unique_ptr<LocalSetI> to_local(const Tensor &d) const override;
    using MultipoleSetI::operator=;
};

extern template class MultipoleSet<2>;

