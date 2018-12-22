// Automatically translated using m2cpp 2.0 on 2018-12-13 19:10:32



#include "LoadGaussFile.h"

/*struct _Property
{
  TYPE Type, degree ;
} ;

TYPE LoadGaussFile(_Property Property)
{
  TYPE GaussDegree, file ;
  _Property Property ;
  if ((strcmp(Property.Type, "1D")))
  {
    GaussDegree = 2*str2num(Property.degree) ;
    file = strcat("Data2/n", num2str(GaussDegree)) ;
  }
  else if ((strcmp(Property.Type, "Triangle")))
  {
    TYPE _var_TYPE = (Property.degree) ;
    if ("1" == _var_TYPE)
    {
      GaussDegree = 1 ;
    }
    else if ("2" == _var_TYPE)
    {
      GaussDegree = 3 ;
    }
    else if ("3" == _var_TYPE)
    {
      GaussDegree = 4 ;
    }
    file = strcat("DataTriangle2/n", num2str(GaussDegree)) ;
  }
  else if ((strcmp(Property.Type, "Quadrilateral")))
  {
    TYPE _var_TYPE = (Property.degree) ;
    if ("1" == _var_TYPE)
    {
      GaussDegree = 2 ;
    }
    else if ("Biquadratic" == _var_TYPE)
    {
      GaussDegree = 3 ;
    }
    else if ("Serendipity" == _var_TYPE)
    {
      GaussDegree = 3 ;
    }
    file = strcat("Data2/n", num2str(GaussDegree)) ;
  }
  else if ((strcmp(Property.Type, "Tetrahedral")))
  {
    TYPE _var_TYPE = (Property.degree) ;
    if ("1" == _var_TYPE)
    {
      GaussDegree = 1 ;
    }
    else if ("2" == _var_TYPE)
    {
      GaussDegree = 4 ;
    }
    else if ("3" == _var_TYPE)
    {
      GaussDegree = 5 ;
    }
    file = strcat("DataTetrahedral2/n", num2str(GaussDegree)) ;
  }
  else if ((strcmp(Property.Type, "Hexahedral")))
  {
    TYPE _var_TYPE = (Property.degree) ;
    if ("1" == _var_TYPE)
    {
      GaussDegree = 2 ;
    }
    else if ("2" == _var_TYPE)
    {
      GaussDegree = 3 ;
    }
    else if ("3" == _var_TYPE)
    {
      GaussDegree = 5 ;
    }
    file = strcat("Data2/n", num2str(GaussDegree)) ;
  }
  return file ;
}*/
std::string FEMtools::LoadGaussFile(const libGmshReader::MeshReader &Mesh, const int &ElementType)
{
    //std::string &ElementName =Mesh.GmshElementName[ElementType];
    int PosOfSpace= Mesh.GmshElementName[ElementType].find_first_of(' ');
    std::string ElementName=Mesh.GmshElementName[ElementType].substr(0,PosOfSpace);
    std::cout<<"Loading Gauss File of Element Type: "<<ElementName<<"\n";
    int &MeshType= Mesh.GmshElementType[ElementType];
    int &p= Mesh.order[ElementType];
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
