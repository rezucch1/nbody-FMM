#include "leaf.hpp"
#include "multipole_set.hpp"
Leaf::Leaf(const Node *parent, const std::vector<Particle>::iterator &particles_begin,
    const std::vector<Particle>::iterator &particles_end): NodeI(parent), particles_begin(particles_begin), particles_end(particles_end)
    {
        calculateMC();
        unsigned int L=2; //to be defined
        multipole_set = std::make_unique<MultipoleSet<2>>(L);
        for( auto p = particles_begin; p< particles_end; ++p){
            Tensor r = p->get_position() - mass_center;
            
            std::unique_ptr<PowerSetI> z;
            if (r.dim == 2)
                z = std::make_unique<PowerSet<2>>(L, r);
            else if (r.dim == 3){
                z = std::make_unique<PowerSet<3>>(L, r);
            }

            (*z) *= p->get_mass();
            (*multipole_set) += *z;
        }


    }

    

    void Leaf::calculateMC()
    {
        mass_center= 0 * particles_begin->get_position();
        double total_mass = 0;
        for(std::vector<Particle>::iterator p = particles_begin; p < particles_end; p++){
            total_mass += p->get_mass();
            mass_center = p->get_mass() * p->get_position();
        }
        mass_center /= total_mass;

        return;
    }
