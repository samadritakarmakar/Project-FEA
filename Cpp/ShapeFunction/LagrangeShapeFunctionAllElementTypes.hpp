#ifndef SHAPEFUNCTIONALLELEMENTTYPES_HPP
#define SHAPEFUNCTIONALLELEMENTTYPES_HPP
#include <armadillo>
#include <vector>
#include "libGmshReader.h"
#include "LagrangeShapeFunction.hpp"
using namespace arma;
class LagrangeShapeFunctionAllElementTypes
{
public:
    /// Constructor
    LagrangeShapeFunctionAllElementTypes(libGmshReader::MeshReader &Mesh)
    {
        //std::vector<ShapeFunction> N(Mesh.NumOfElementTypes);
        ShpFnc=std::vector<LagrangeShapeFunction>(Mesh.NumOfElementTypes);
        N=std::vector<mat>(Mesh.NumOfElementTypes);
        for (int i=0;i<Mesh.NumOfElementTypes;i++)
        {
            ShpFnc[i]=LagrangeShapeFunction(Mesh,i);
        }
    }

     std::vector<mat> GetShapeFunction(mat &GaussPointx)
    {
        for (int i=0; i<ShpFnc.size(); i++)
        {
            N[i]=ShpFnc[i].LagrangeShapeFunction::GetShapeFunction(GaussPointx);
        }
        return N;
    }
     std::vector<mat> GetShapeFunction(mat &GaussPointx, mat &GaussPointy)
    {
        for (int i=0; i<ShpFnc.size(); i++)
        {
            N[i]=ShpFnc[i].LagrangeShapeFunction::GetShapeFunction(GaussPointx, GaussPointy);
        }
        return N;
    }
     std::vector<mat> GetShapeFunction(mat &GaussPointx, mat &GaussPointy, mat &GaussPointz)
    {
        for (int i=0; i<ShpFnc.size(); i++)
        {
            N[i]=ShpFnc[i].LagrangeShapeFunction::GetShapeFunction(GaussPointx, GaussPointy, GaussPointz);
        }
        return N;
    }

protected:
std::vector<LagrangeShapeFunction> ShpFnc;
std::vector<mat> N;

};
#endif // SHAPEFUNCTIONALLELEMENTTYPES_HPP
