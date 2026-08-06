#pragma once
#include <tuple>

std::tuple<double, double> mult_complex(double Re1, double Im1, double Re2, double Im2){
    return { Re1*Re2 - Im1*Im2, Re1*Im2 + Re2*Im1};
        
}
