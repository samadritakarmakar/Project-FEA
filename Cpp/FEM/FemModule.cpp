#include "FemModule.h"
void FemModule::SetNoOfInstances(const libGmshReader::MeshReader &Mesh, const Models &model)
{
    GaussLength=new int [Mesh.NumOfElementTypes];
   DependentEpsilon=new int [Mesh.NumOfElementTypes];
   //dof=new int [Mesh.NumOfElementTypes];
   NodePositions=new umat [model.coupleLevel+1];
}
FemModule::~FemModule()
{
    delete []GaussLength;
    delete []DependentEpsilon;
    delete []dof;
    delete []NodePositions;
}

void FemModule::SetNodePositions(int couplingNumber, const umat &ElementNodes, int vectorLevel)
{
    GetNodePostions(NodePositions[couplingNumber], ElementNodes, vectorLevel);
}

FemModule::FemModule(const libGmshReader::MeshReader &Mesh,  const Models &model)
{
    SetNoOfInstances(Mesh, model);
    const int &maxNodeNumber= Mesh.maxNodeNumber;
    const uvec &VectorLevel=model.VectorLevel;
    const int &coupleLevel=model.coupleLevel;
    int dofFinalGlobal=maxNodeNumber*accu(VectorLevel); //To be used to set the size of Final LHS Global and Final RHS Global Matrices
    uvec dofPerTerm=maxNodeNumber*(VectorLevel);
    int NoOfMatrixTerms=(coupleLevel+1)*(coupleLevel+1); //To be used to set number of LHS and RHS matrices
    for (int i=0; i<Mesh.NumOfElementTypes; i++)
    {
        int NumOfElements=Mesh.ElementNodes[i].n_rows;
        std::string GaussFileName =FEMtools::LoadGaussFile(Mesh, i);
        mat data;
        data.load(GaussFileName);
        GaussLength[i]=data.n_rows; //Number of Gauss Points needed per intergration
        DependentEpsilon[i]=data.n_cols;
        DependentEpsilon[i]=DependentEpsilon[i]-2; //Number of Gauss points per weight. Factoring for the numbering and the weights in the Gauss Files.
        for (int couplingNumber=0; couplingNumber<coupleLevel+1; couplingNumber++)
        {
            SetNodePositions(couplingNumber, Mesh.ElementNodes[i], VectorLevel(couplingNumber));
            //std::cout<<"Max Node Pos= "<<max(vectorise(NodePositions[couplingNumber]))<<"\n";
            //std::cout<<"Min Node Pos= "<<min(vectorise(NodePositions[couplingNumber]))<<"\n";
        }
    }
}


