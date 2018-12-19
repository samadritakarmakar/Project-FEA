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
std::string FEMtools::LoadGaussFile(const libGmshReader::MeshReader &Mesh, const int & ElementType)
{
    int &MeshType= Mesh.GmshElementType[ElementType];
    int &p= Mesh.order[ElementType];
    float n;

    std::string GaussFileName="GaussData/";
    //for line elements
    if (MeshType == 1 || MeshType == 8 || MeshType == 26 || MeshType == 27 || MeshType == 28 )
    {
        FEMtools::TensorGauss(p, n, GaussFileName);
    }
    //for triangle elements
    else if (MeshType == 2)
    {
        n=1;
        GaussFileName=GaussFileName+"DataTriangle/n"+std::to_string(int (n));
    }
    else if (MeshType == 9)
    {
        n=3;
        GaussFileName=GaussFileName+"DataTriangle/n"+std::to_string(int (n));
    }
    else if (MeshType == 20  || MeshType == 21 || MeshType == 22 || MeshType == 23)
    {
        n=6;
        GaussFileName=GaussFileName+"DataTriangle/n"+std::to_string(int (n));
    }
    else if (MeshType == 24 || MeshType == 25)
    {
        n=7;
        GaussFileName=GaussFileName+"DataTriangle/n"+std::to_string(int (n));
    }
    //For Quadrilaterals
    else if (MeshType == 3 || MeshType == 10 || MeshType == 16)
    {
        FEMtools::TensorGauss(p, n, GaussFileName);
    }
    return GaussFileName;
    //
}

void FEMtools::TensorGauss(int &p, float &n, std::string &GaussFileName)
{
    n=(p+2.0)/(float (p));
    n=std::ceil(n);
    GaussFileName=GaussFileName+"Data/n"+std::to_string(int (n));
}
