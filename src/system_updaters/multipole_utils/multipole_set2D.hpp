#pragma once

#include "multipole_set.hpp"

template<>
class MultipoleSet<2> : public MultipoleSetI {

    protected:
    virtual std::ostream &output_to_stream(std::ostream &os) const override{
        os << (*this)(0);
        for (uint l = 1; l <= L; ++l)
            os << ", " << (*this)(l);
        return os;
    }
    
    public: 
    MultipoleSet(unsigned int);
    virtual std::complex<double> operator()(unsigned l, int m=0) const override;
    virtual std::unique_ptr<MultipoleSetI> weigh_children_with_distance(const Tensor &d) const override;
};

extern template class MultipoleSet<2>;

