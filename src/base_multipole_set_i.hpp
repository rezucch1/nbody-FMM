#pragma once
#include <vector>
#include <complex>



class BaseMultipoleSetI{
    protected :
    std::vector<double> elements;
    unsigned int L;
    
    BaseMultipoleSetI(unsigned int L);
    public :
    virtual std::complex<double> operator()(unsigned l, int m=0) const=0;
    
    friend class MultipoleSetI;

} ;