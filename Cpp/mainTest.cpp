#include<iostream>
#include "libMeshReader/libMeshReader.h"
#include <armadillo>
using namespace arma;
int main(int argc, char *argv[])
{
    MeshReader a;
    mat sortedVariable, unsortedVariable;
    unsortedVariable={{1,36,43},{36,65,43}};
    sortedVariable= a.SortAndKeep(unsortedVariable);
    std::cout<<sortedVariable;
    return 0;
}
