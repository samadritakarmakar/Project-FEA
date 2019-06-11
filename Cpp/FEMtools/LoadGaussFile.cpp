// Automatically translated using m2cpp 2.0 on 2018-12-13 19:10:32


#include "LoadGaussFile.h"

std::string FEMtools::LoadGaussFile(const libGmshReader::MeshReader &Mesh, const int &ElementType)
{
    //std::string &ElementName =Mesh.GmshElementName[ElementType];
    int PosOfSpace= Mesh.GmshElementName[ElementType].find_first_of(' ');
    std::string ElementName=Mesh.GmshElementName[ElementType].substr(0,PosOfSpace);
    std::cout<<"Loading Gauss File of Element Type: "<<ElementName<<"\n";
    int MeshType= Mesh.GmshElementType[ElementType];
    int p= Mesh.order[ElementType];
    float n;

    std::string GaussFileName="GaussData/";
    //for line elements
    if (ElementName.compare("Line")==0)
    {
        FEMtools::TensorGauss(p, n, GaussFileName);
    }
    //for triangle elements
    else if (ElementName.compare("Triangle")==0)
    {
        if (p==1)
        {
            n=1;
        }
        else if (p==2)
        {
            n=3;
        }
        else if (p==3 || p==4)
        {
            n=6;
        }
        else
        {
            n=7;
        }
        GaussFileName=GaussFileName+"DataTriangle/n"+std::to_string(int (n));
    }
    //For Quadrilaterals
    else if (ElementName.compare("Quadrilateral")==0)
    {
        FEMtools::TensorGauss(p, n, GaussFileName);
    }
    //For Tetrahedrals
    else if (ElementName.compare("Tetrahedron")==0)
    {
        if (p==1)
        {
            n=1;
        }
        else if (p==2)
        {
            n=4;
        }
        else {
            n=15;
        }
        GaussFileName=GaussFileName+"DataTetrahedral/n"+std::to_string(int (n));
    }
    //For Hexahedrals
    else if (ElementName.compare("Hexahedron")==0)
    {
        FEMtools::TensorGauss(p, n, GaussFileName);
    }
    else if (ElementName.compare("Prism")==0)
    {
        FEMtools::TensorGauss(p, n, GaussFileName);
    }
    else if (ElementName.compare("Pyramid")==0)
    {
        FEMtools::TensorGauss(p, n, GaussFileName);
    }
    else
    {
        std::cout<<"Could not find appropriate Quadrature\n";
    }
    std::cout<<"Loading Quadrature File: "<<GaussFileName<<"\n";
    return GaussFileName;
    //
}

void FEMtools::TensorGauss(int &p, float &n, std::string &GaussFileName)
{
    n=(p+2.0)/(2.0);
    n=std::ceil(n);
    GaussFileName=GaussFileName+"Data/n"+std::to_string(int (n));
}
