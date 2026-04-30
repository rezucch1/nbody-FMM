#include "multipole.hpp"
constexpr unsigned int threshold = 5;

Multipole::Multipole(const Multipole *parent, std::vector<Particle> particles, const Tensor &a, const Tensor &b)
: parent(parent)
{
    const Tensor c = (a+b)/2;
    // unsigned int counter[2][2]; //2D

    // for(int i=0; i<2; i++){
    //     for(int j=0; j<2; j++){
    //         counter[i][j]=0;
    //     }
    // }

    std::vector<Particle> childrenparts[2][2];
    int i, j;
    for(Particle &p: particles){
        i = p.get_position()[0] > c[0]; //if TRUE, i=1; if FALSE, i=0;
        j = p.get_position()[1] > c[1]; //if TRUE, j=1; if FALSE, j=0;

        // counter[i][j]++;
        childrenparts[i][j].push_back(p);
    }


    for( i=0; i<2; i++){
        for( j=0; j<2; j++){
            Tensor a2(a.dim), b2(b.dim);

            if( childrenparts[i][j].size()> threshold) //non è leaf
            {
                if(i==0){
                    a2[1]=a[1];
                    b2[1]=c[1];
                }
                else{
                    a2[1]=c[1];
                    b2[1]=b[1];
                }

                if(j==0){
                    a2[0]=a[0];
                    b2[0]=c[0];
                }
                else{
                    a2[0]=c[0];
                    b2[0]=b[0];
                }

                children.push_back(Multipole(this, childrenparts[i][j], a2, b2));
            }
            else //è leaf
            {

            }
        }
    }

}