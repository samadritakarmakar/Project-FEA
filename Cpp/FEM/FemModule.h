#ifndef FEMMODULE_H
#define FEMMODULE_H
#include <armadillo>
#include "libGmshReader.h"
#include "LoadGaussFile.h"
#include "models.h"
using namespace arma;
class FemModule
{
public:

    FemModule(libGmshReader::MeshReader &Mesh, Models &model);
    int *GaussLength,  *DependentEpsilon, *dof;
    ~FemModule();
protected:
    void SetNoOfInstances(libGmshReader::MeshReader &Mesh);
};

#endif // FEMMODULE_H
