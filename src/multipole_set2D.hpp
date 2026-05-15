#pragma once

#include "multipole_set.hpp"

template<>
class MultipoleSet<2> : public MultipoleSetI {
    
    public: 
    MultipoleSet(unsigned int L, const Tensor &r);

};

