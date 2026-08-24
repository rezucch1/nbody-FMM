#pragma once

#include "power_set.hpp"

class Irregular : public PowerSetI {
    
    public: 
    Irregular(unsigned int L, const Tensor &r);
    virtual std::complex<double> operator()(unsigned l, int m=0) const override;

};

