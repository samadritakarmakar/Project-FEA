#include "FemModule.h"

FemModule::FemModule(libGmshReader::MeshReader &Mesh)
{
    for (int i=0; i<Mesh.NumOfElementTypes; i++)
    {
        int NumOfElements=Mesh.ElementNodes[i].n_rows;
        std::string GaussFileName =FEMtools::LoadGaussFile(Mesh, i);
        mat data;
        data.load(GaussFileName);
    }
}
