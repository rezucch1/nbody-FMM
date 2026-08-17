#include "leaf.hpp"
#include "multipole_set.hpp"
Leaf::Leaf(const Node *parent, Particle** particles_begin,
    Particle** particles_end): NodeI(parent), particles_begin(particles_begin), particles_end(particles_end)
    {
        calculateMC();
        unsigned int L=2; //to be defined
        multipole_set = std::make_unique<MultipoleSet<2>>(L);
        (*multipole_set) = 0;
        for( auto p = particles_begin; p< particles_end; ++p){
            Tensor r = (*p)->get_position() - mass_center;
            
            std::unique_ptr<PowerSetI> z;
            if (r.dim == 2)
                z = std::unique_ptr<PowerSet<2>>(new PowerSet<2>(L, r));
            // else if (r.dim == 3){
            //     z = std::make_unique<PowerSet<3>>(L, r);
            // }

            (*z) *= (*p)->get_mass();
            (*multipole_set) += *z;
        }


    }

    

    void Leaf::calculateMC()
    {
        double total_mass = 0;

        // The mass_center should be initialized anyways, at least to 0;
        mass_center = 0 * (*particles_begin)->get_position();

        if (particles_begin == particles_end)
            return;

        for(Particle** p = particles_begin; p < particles_end; p++){
            total_mass += (*p)->get_mass();
            mass_center += (*p)->get_mass() * (*p)->get_position();
        }
        mass_center /= total_mass;

        return;
    }
