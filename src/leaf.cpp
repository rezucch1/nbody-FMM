#include "leaf.hpp"
#include "multipole_set.hpp"
Leaf::Leaf(const Node *parent, const std::vector<Particle>::iterator &particles_begin,
    const std::vector<Particle>::iterator &particles_end): Multi_inter(parent), particles_begin(particles_begin), particles_end(particles_end)
    {
        calculateMC();
        unsigned int L=2; //to be defined
        std::unique_ptr<MultipoleSetI> M = std::make_unique<MultipoleSet<2>>(L);
        for( auto p = particles_begin; p< particles_end; ++p){
            Tensor r= p->get_position() - mass_center;
            
            std::unique_ptr<PowerSetI> z;
            if (r.dim == 2)
             z = std::make_unique<PowerSet<2>>(L, r);
            else if (r.dim == 3){
                //TODO
            }

            (*z) *= p->get_mass();
            (*M) += *z;
        }


    }

    

    void Leaf::calculateMC()
    {
        mass_center={0.0, 0.0};
        double total_mass = 0;
        //2D
        for(std::vector<Particle>::iterator p = particles_begin; p < particles_end; p++){
            total_mass += p->get_mass();
            mass_center = p->get_mass() * p->get_position();
        }
        mass_center /= total_mass;

        return;
    }
