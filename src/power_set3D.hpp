#pragma once

#include "power_set.hpp"

template<>
class PowerSet<3> : public PowerSetI {
    
    public: 
    PowerSet(unsigned int L, const Tensor &r);
    virtual std::complex<double> operator()(unsigned l, int m) override;

};

