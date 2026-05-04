#include "multipole.hpp"
constexpr unsigned int threshold = 5;

Multipole::Multipole(const Multipole *parent, std::vector<Particle>::iterator &particles_begin, std::vector<Particle>::iterator &particles_end, const Tensor &a, const Tensor &b)
: parent(parent)
{
    const Tensor c = (a+b)/2;
    unsigned int counter[4]; //2D

    for(int i=0; i<4; i++){
            counter[i]=0; //initialize all counter to 0
    }

    std::vector<Particle>::iterator children_id[5];
    int i;
    for(std::vector<Particle>::iterator p = particles_begin; p < particles_end; p++){
        i = p->get_position()[0] > c[0]; //c[0] is the x position of c
        //if particle p is on c's left, i=0; if particle p is on c's right, i=1;
        i += 2 * (p->get_position()[1] > c[1]); //if TRUE, i+=10; if FALSE i=i;
        //we add 10 if particle p is below c (in y terms)
        counter[i]++; //count how many particles inside the child multiltpole

    }

    children_id[0] = particles_begin;

    for(i = 1; i < 5 ; i++){
        children_id[i] = children_id[i - 1] + counter[i - 1];
        //each cell contains the address to the first particle of each child
        //but we still need to reorder the array
    }

    //let's STACK
    {
        //we create a temporary array pf particles
        Particle part_temp[particles_end - particles_begin];
        std::move(particles_begin, particles_end, part_temp); //move all our particles array in temporary array
        //now our initial particles array is EMPTY

        for(i=0; i<4; i++){
            counter[i]=0; //initialize all counter to 0
        }
        for(std::vector<Particle>::iterator p = particles_begin; p < particles_end; p++){
            i = p->get_position()[0] > c[0]; //c[0] is the x position of c
            //if particle p is on c's left, i=0; if particle p is on c's right, i=1;
            i += 2 * (p->get_position()[1] > c[1]); //if TRUE, i+=10; if FALSE i=i;
            //we add 10 (binary, so it's 2) if particle p is below c (in y terms)

            *(children_id[i] + counter[i]) = std::move(*p);
            //children_id refers to the start of the child-cell i
            //counter[i] refers to the offset
            counter[i]++; //count how many particles inside the child multiltpole
        }

    }


    for( i=0; i<4; i++){
        
        Tensor a2(a.dim), b2(b.dim);

        if( counter[i] > threshold) //non è leaf
        {
            if(i & 1 == 0){//if the least significant bit=0
                a2[1]=a[1];
                b2[1]=c[1];
            }
            else{
                a2[1]=c[1];
                b2[1]=b[1];
            }

            if(i & 2 == 0){ //if the second least significant bit=0
                a2[0]=a[0];
                b2[0]=c[0];
            }
            else{
                a2[0]=c[0];
                b2[0]=b[0];
            }

            children.push_back(Multipole(this, children_id[i], children_id[i + 1] , a2, b2));
        }
        else //è leaf
        {

        }
        
    }

}