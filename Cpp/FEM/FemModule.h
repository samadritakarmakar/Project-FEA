#ifndef FEMMODULE_H
#define FEMMODULE_H
#include <armadillo>
#include "libGmshReader.h"
//#include "LoadGaussFile.h"
#include "FEMtools.h"
#include "models.h"
using namespace arma;
class FemModule
{
public:

    FemModule(const libGmshReader::MeshReader &Mesh, const Models &model);
    int *GaussLength=nullptr,  *DependentEpsilon=nullptr, *dof=nullptr;
    umat *NodePositions=nullptr;
    ~FemModule();
protected:
    void SetNoOfInstances(const libGmshReader::MeshReader &Mesh, const Models &model);
    void SetNodePositions(int couplingNumber, const umat &ElementNodes, int vectorLevel);
    void AssemblePerTermMatrix();
};

#endif // FEMMODULE_H
