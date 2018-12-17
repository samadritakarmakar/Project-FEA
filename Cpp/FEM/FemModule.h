#ifndef FEMMODULE_H
#define FEMMODULE_H
#include <armadillo>
#include "libGmshReader.h"

using namespace arma;
class FemModule
{
public:
    FemModule(libGmshReader::MeshReader Mesh);
};

#endif // FEMMODULE_H
