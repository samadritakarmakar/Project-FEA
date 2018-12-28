#include "FemModule.h"
void FemModule::SetNoOfInstances(libGmshReader::MeshReader &Mesh)
{
    GaussLength=new int [Mesh.NumOfElementTypes];
   DependentEpsilon=new int [Mesh.NumOfElementTypes];
   dof=new int [Mesh.NumOfElementTypes];
}
FemModule::~FemModule()
{
    delete []GaussLength;
    delete []DependentEpsilon;
    delete []dof;
}

FemModule::FemModule(libGmshReader::MeshReader &Mesh, Models &model)
{
    SetNoOfInstances(Mesh);
    for (int i=0; i<Mesh.NumOfElementTypes; i++)
    {
        int NumOfElements=Mesh.ElementNodes[i].n_rows;
        std::string GaussFileName =FEMtools::LoadGaussFile(Mesh, i);
        mat data;
        data.load(GaussFileName);
        GaussLength[i]=data.n_rows;
        DependentEpsilon[i]=data.n_cols;
        DependentEpsilon[i]=DependentEpsilon[i]-2; //Factoring for the numbering and the weights in the Gauss Files
        dof[i]=Mesh.ContainsNodes[i].n_rows*model.VectorLevel;
    }
}
