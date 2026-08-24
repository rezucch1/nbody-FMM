#pragma once

#include "multipole_set.hpp"

template<>
class MultipoleSet<3> : public MultipoleSetI {
    
    public: 
    MultipoleSet(unsigned int L);
    virtual std::complex<double> operator()(unsigned l, int m) const override;
    virtual MultipoleSetI &operator+=(const PowerSetI &z) override;
    virtual std::unique_ptr<MultipoleSetI> weigh_children_with_distance (const Tensor &d) const override;
    virtual std::unique_ptr<LocalSetI> to_local(const Tensor &d) const;
};

