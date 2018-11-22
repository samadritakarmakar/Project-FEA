#include<iostream>
#include "libMeshReader/libMeshReader.h"
#include <armadillo>
using namespace arma;
int main()
{
    MeshReader a;
    mat sortedVariable, unsortedVariable;
    unsortedVariable={{1,36,43},{36,65,43}};
    sortedVariable= a.SortAndKeep(unsortedVariable);
    std::cout<<sortedVariable;

    mat b,c,d;
    d=zeros(3,1);
    b={{1,2,3},{5,6,7},{7,10,9}};
    c<<1<<endr
    <<2<<endr
    <<3<<endr;
    d=solve(b,c);
    std::cout<<d;
    return 0;



}
