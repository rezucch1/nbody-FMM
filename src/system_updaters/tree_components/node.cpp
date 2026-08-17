#include "node.hpp"
#include "leaf.hpp"
constexpr unsigned int threshold = 5;

Node::Node(const Node *parent, Particle** particles_begin, Particle** particles_end, const Tensor &a, const Tensor &b)
    : NodeI(parent)
{
    const Tensor c = (a + b) / 2;
    unsigned int counter[4]; // 2D, dim

    for(int i=0; i<4; i++){
        counter[i]=0; //initialize all counter to 0
    }

    Particle** children_id[5]; //2^dim + 1 divisione sottovettori punti estremi dei sottovettori
    int i;
    for(Particle** p = particles_begin; p < particles_end; p++){
        i = (*p)->get_position()[0] > c[0]; //c[0] is the x position of c
        //if particle p is on c's left, i=0; if particle p is on c's right, i=1;
        i += 2 * ((*p)->get_position()[1] > c[1]); //if TRUE, i+=10; if FALSE i=i;
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
        Particle* part_temp[particles_end - particles_begin];
        std::move(particles_begin, particles_end, part_temp); //move all our particles array in temporary array
        //now our initial particles array is EMPTY

        for(i=0; i<4; i++){
            counter[i]=0; //initialize all counter to 0
        }

        for(int p = 0; p < particles_end - particles_begin; p++){
            i = part_temp[p]->get_position()[0] > c[0]; //c[0] is the x position of c
            //if particle p is on c's left, i=0; if particle p is on c's right, i=1;
            i += 2 * (part_temp[p]->get_position()[1] > c[1]); //if TRUE, i+=10; if FALSE i=i;
            //we add 10 (binary, so it's 2) if particle p is below c (in y terms)

            *(children_id[i] + counter[i]) = part_temp[p];
            //children_id refers to the start of the child-cell i
            //counter[i] refers to the offset
            counter[i]++; //count how many particles inside the child multiltpole
        }

    }


    for( i=0; i<4; i++){
        
        Tensor a2(a.dim), b2(b.dim);
        std::unique_ptr<NodeI> child;

        if( counter[i] > threshold) //non è leaf
        {
            if((i & 1) == 0){//if the least significant bit=0
                a2[0]=a[0];
                b2[0]=c[0];
            }
            else{
                a2[0]=c[0];
                b2[0]=b[0];
            }

            if((i & 2) == 0){ //if the second least significant bit=0
                a2[1]=a[1];
                b2[1]=c[1];
            }
            else{
                a2[1]=c[1];
                b2[1]=b[1];
            }

            child = std::make_unique<Node>(this, children_id[i], children_id[i + 1] , a2, b2);
        }
        else // è leaf
        {
            child = std::make_unique<Leaf>(this, children_id[i], children_id[i + 1]);
        }

        children.push_back(std::move(child));
    }

    calculateMC();
    unsigned int L=2; //to be defined
    multipole_set = std::make_unique<MultipoleSet<2>>(L);
    (*multipole_set) = 0;
    for (const auto &child : children)
        *multipole_set += *getMultipoleSet(*child)->weigh_children_with_distance(NodeI::getMassCenter(*child) - mass_center);

}

void Node::calculateMC()
{
    mass_center= 0 * getMassCenter(*children[0]);
    double total_mass = 0;
    for(const auto &c : children){
        double child_mass = (*getMultipoleSet(*c))(0).real();
        total_mass += child_mass;
        mass_center = child_mass * NodeI::getMassCenter(*c);
    }
    mass_center /= total_mass;

    return;
}
