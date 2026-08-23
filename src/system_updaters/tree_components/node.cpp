#include "node.hpp"
#include "leaf.hpp"
constexpr unsigned int threshold = 5;

Node::Node(std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int id_child, Particle** particles_begin, Particle** particles_end, const Tensor &a, const Tensor &b)
    : NodeI(allocator, depth, id_child ,particles_begin, particles_end)
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

    if (allocator.size() == depth + 1)
        allocator.emplace_back(1<<(dim * (depth + 1)));


    for( i=0; i < 1<<dim; i++){
        
        Tensor a2(a.dim), b2(b.dim);
        std::unique_ptr<NodeI> child;
        unsigned int current_child_id = (id_child << dim) + i;

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

            child = std::make_unique<Node>(allocator, depth + 1, current_child_id, children_id[i], children_id[i + 1] , a2, b2);
        }
        else if (counter[i] > 0) // è leaf
        {
            child = std::make_unique<Leaf>(allocator, depth + 1, current_child_id, children_id[i], children_id[i + 1]);
        }
        allocator[depth + 1][current_child_id] = std::move(child);
    }
}

// const std::vector<std::unique_ptr<NodeI>> &Node::get_children() const{
//   return children;
// }

void Node::compute_multipoles(unsigned int L){

    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        (*c)->compute_multipoles(L);

    calculateMC();

    multipole_set = std::make_unique<MultipoleSet<2>>(L);
    (*multipole_set) = 0;
    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        *multipole_set += *getMultipoleSet(**c)->weigh_children_with_distance(NodeI::getMassCenter(**c) - mass_center);
}

void Node::collect_multipoles_to_locals(){
    NodeI::collect_multipoles_to_locals();

    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        (*c)->collect_multipoles_to_locals();
}

void Node::propagate_locals(LocalSetI *parent_local){
    NodeI::propagate_locals(parent_local);

    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        (*c)->propagate_locals(local_set->distribute_parent_with_distance(getMassCenter(**c) - mass_center));
}

const std::vector<NodeI *> &Node::get_neighbours() const
{
  return neighbours_list;
}

void Node::get_partition(std::vector<std::tuple<Particle *, int, int>> &partitions) const{
    NodeI::get_partition(partitions);

    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        (*c)->get_partition(partitions);
}

NodeI *Node::get_child(unsigned int idx)
{
  return allocator[depth + 1][id_child * 1<<dim + idx].get();
}

std::vector<std::unique_ptr<NodeI>>::iterator Node::get_children_begin(){
  return allocator[depth + 1].begin() + (id_child << dim);
}

std::vector<std::unique_ptr<NodeI>>::iterator Node::get_children_end(){
  return allocator[depth + 1].begin() + ((id_child+1) << dim);
}

std::vector<std::unique_ptr<NodeI>>::const_iterator Node::get_children_begin() const
{
  return allocator[depth + 1].cbegin() + (id_child << dim);
}

std::vector<std::unique_ptr<NodeI>>::const_iterator Node::get_children_end() const
{
  return allocator[depth + 1].cbegin() + ((id_child+1) << dim);
}

void Node::calculateMC()
{
    mass_center = Tensor(dim);
    mass_center *= 0;
    double total_mass = 0;
    for(auto c = get_children_begin(); c < get_children_end(); ++c) if (*c){
        double child_mass = (*getMultipoleSet(**c))(0).real();
        total_mass += child_mass;
        mass_center = child_mass * NodeI::getMassCenter(**c);
    }
    mass_center /= total_mass;

    return;
}
