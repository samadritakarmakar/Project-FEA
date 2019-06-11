#ifndef SHAPEFUNCTION_H
#define SHAPEFUNCTION_H
#include "libGmshReader.h"
#include <armadillo>
#include <string>
#include "FindCoeffLine.hpp"
#include "FindCoeffTriangle.hpp"
#include "FindCoeffQuad.hpp"
#include "FindCoeffTertrahedral.hpp"
#include "FindCoeffHexahedral.hpp"
using namespace arma;
class LagrangeShapeFunction
{
public:
    mat (LagrangeShapeFunction::*FunctionPtr)( mat&, mat&, mat&);
    /// Default Constructor. Remember to call 'ShapeFunction(libGmshReader::MeshReader &Mesh, const int& ElementType)'
    LagrangeShapeFunction()
    {
        //----Do Nothing-----
    }
    /// Constructor. Here type of element is checked and the needed Function pointer is assigned.
    /// The Shape Functions for each Gauss Point is row wise.
    /// [N1G1, N1G2, N1G3,...; N2G1, N2G2, N2G3,...; N3G1, N3G2, N3G3,...;...]
    LagrangeShapeFunction(libGmshReader::MeshReader &Mesh, const int& ElementType)
    {
        int i=ElementType;
        order=Mesh.order[ElementType];
        if(Mesh.GmshElementNameOnly[i].compare("Line")==0)
        {
            FunctionPtr=&LagrangeShapeFunction::LineShapeFuntion;
            Coefficient=FindCoeffLine(order).t();
        }
        else if(Mesh.GmshElementNameOnly[i].compare("Triangle")==0)
        {
            FunctionPtr=&LagrangeShapeFunction::TriangleShapeFunction;
            Coefficient=FindCoeffTriangle(order).t();
        }
        else if(Mesh.GmshElementNameOnly[i].compare("Quadrilateral")==0)
        {
            FunctionPtr=&LagrangeShapeFunction::QuadShapeFunction;
            Coefficient=FindCoeffQuad(order).t();
        }
        else if(Mesh.GmshElementNameOnly[i].compare("Tetrahedron")==0)
        {
            FunctionPtr=&LagrangeShapeFunction::TetrahedralShapeFunction;
            Coefficient=FindCoeffTetrahedral(order).t();
        }
        else if(Mesh.GmshElementNameOnly[i].compare("Hexahedron")==0)
        {
            FunctionPtr=&LagrangeShapeFunction::HexahedralShapeFunction;
            Coefficient=FindCoeffHexadedral(order).t();
        }
        else
        {
            std::cout<<"Shape Function "<<Mesh.GmshElementNameOnly[i]<<" is not supported!!\n";
            throw;
        }
    }
    ///For 1D Shape Functions
    /// The Shape Functions for each Gauss Point is row wise.
    /// [N1G1, N1G2, N1G3,...; N2G1, N2G2, N2G3,...; N3G1, N3G2, N3G3,...;...]
    inline mat GetShapeFunction( mat &GaussPointx)
    {
        mat GaussPointz={0}; mat GaussPointy={0};
        return (this->*FunctionPtr)( GaussPointx, GaussPointy, GaussPointz);
    }
    ///For 2D Shape Functions
    /// The Shape Functions for each Gauss Point is row wise.
    /// [N1G1, N1G2, N1G3,...; N2G1, N2G2, N2G3,...; N3G1, N3G2, N3G3,...;...]
    inline mat GetShapeFunction( mat &GaussPointx, mat &GaussPointy)
    {
        mat GaussPointz={0};
        return (this->*FunctionPtr)( GaussPointx, GaussPointy, GaussPointz);
    }
    ///For 3D Shape Functions
    /// The Shape Functions for each Gauss Point is row wise.
    /// [N1G1, N1G2, N1G3,...; N2G1, N2G2, N2G3,...; N3G1, N3G2, N3G3,...;...]
    inline mat GetShapeFunction( mat &GaussPointx, mat &GaussPointy, mat &GaussPointz)
    {
        return (this->*FunctionPtr)( GaussPointx, GaussPointy, GaussPointz);
    }


private:
    inline mat LineShapeFuntion( mat &GaussPointx, mat &GaussPointy, mat &GaussPointz)
    {
        mat Matrx;
        BuildLineMatrix(order, GaussPointx, NumOfLineNodes(order), Matrx);
        mat N=(Matrx*Coefficient).t();
        return N;
    }
    inline mat TriangleShapeFunction( mat &GaussPointx, mat &GaussPointy, mat &GaussPointz)
    {
        mat Matrx;
        int NumOfNodes=NumOfTriNodes(order);
        BuildTriangleMatrix(order, GaussPointx, GaussPointy, NumOfNodes, Matrx);
        mat N=(Matrx*Coefficient).t();
        return N;
    }
    inline mat QuadShapeFunction( mat &GaussPointx, mat &GaussPointy, mat &GaussPointz)
    {
        mat Matrx;
        int NumOfNodes=NumOfQuadNodes(order);
        BuildQuadMatrix(order, GaussPointx, GaussPointy, NumOfNodes, Matrx);
        mat N=(Matrx*Coefficient).t();
        return N;
    }

    inline mat TetrahedralShapeFunction( mat &GaussPointx, mat &GaussPointy, mat &GaussPointz)
    {
        mat Matrx;
        int NumOfNodes=NumOfTetraNodes(order);
        BuildTetrahedralMatrix(order, GaussPointx, GaussPointy, GaussPointz, NumOfNodes, Matrx);
        mat N=(Matrx*Coefficient).t();
        return N;
    }

    inline mat HexahedralShapeFunction( mat &GaussPointx, mat &GaussPointy, mat &GaussPointz)
    {
        mat Matrx;
        int NumOfNodes=NumOfHexaNodes(order);
        BuildHexahedralMatrix(order, GaussPointx, GaussPointy, GaussPointz, NumOfNodes, Matrx);
        mat N=(Matrx*Coefficient).t();
        return N;
    }
    mat Coefficient;
    int order;

};

#endif // SHAPEFUNCTION_H
