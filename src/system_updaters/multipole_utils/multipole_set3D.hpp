#pragma once

#include "multipole_set.hpp"

template<>
class MultipoleSet<3> : public MultipoleSetI {
    
    public: 
    MultipoleSet(unsigned int L);
    virtual std::complex<double> operator()(unsigned l, int m) override;

};

