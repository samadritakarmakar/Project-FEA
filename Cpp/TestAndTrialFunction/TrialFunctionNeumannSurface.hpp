#ifndef TRIALFUNCTIONNEUMANNSURFACE_HPP
#define TRIALFUNCTIONNEUMANNSURFACE_HPP

#include <TrialFunction.hpp>
#include "libGmshReader.h"

class TrialFunctionNeumannSurface: public TrialFunction
{
public:
    TrialFunctionNeumannSurface (TrialFunction& u, int PhysicalGroupNumber): TrialFunction (*u.Msh, 2, u.vectorLvl-1),
        originalDimension(u.Msh->ElementData::dim), PhysclGrpNum (PhysicalGroupNumber)
    {
        //cout<<"Grad(u) Neumann Surface\n"<<mat(Get_grad_u(0,0,0));
        //cout<<"dS_norm = "<<Get_dS(0,0,0)<<"\n";
        GetNumberOfVariables();
        originalVctrLvl=(u.vectorLvl);
        Regenerate_u();
    }

    /// Sets the NoOfElements and ElmntNodes or 'Connectivity Matrix' for a certain Element Type.
    void GetNumberOfVariables()
    {
        for (int ElementType = 0; ElementType<Msh->NumOfElementTypes; ++ElementType)
        {
            if(int(Msh->ElmntPhysclGrpNodes[ElementType].size())==0)
            {
                cout<<"A Physical Group at Index "<<PhysclGrpNum<<" does not exist!!!"<<"\n";
                throw;
            }
            NoOfElements[ElementType]= Msh->ElmntPhysclGrpNodes[ElementType][PhysclGrpNum].n_rows;
            ElmntNodes[ElementType]= Msh->ElmntPhysclGrpNodes[ElementType][PhysclGrpNum];
        }
    }
    /// This function generates the F matrix. This in mathematical terms is the Jacobian dx/dEps
    void Get_F(int ElementType, int ElementNumber, int GaussPntr, mat& F)
    {
        umat NodesAtElmntNmbr=Msh->ElmntPhysclGrpNodes[ElementType][PhysclGrpNum].row(ElementNumber);
        //cout<<"Gmsh Node Tags are "<<Msh->GmshNodeTag[ElementType].row(ElementNumber);
        //cout<<"ElementNodes =\n"<<NodesAtElmntNmbr;
        //cout<<"Coodinates of "<<ElementNumber<<" are \n"<<Msh->NodalCoordinates.rows(NodesAtElmntNmbr);
        mat Coordinates=Msh->NodalCoordinates.rows(NodesAtElmntNmbr);
        //cout<<"Coodinates of "<<ElementNumber<<" are \n"<<Coordinates.cols(0,MeshDimension-1)<<"\n";
        //coords of x for dim 1; x & y for dim 2; x, y & z for dim 3;
        F=Coordinates.cols(0,MeshDimension-1).t()*dN_by_dEps[ElementType][GaussPntr];
        //cout<<"Coordinates.cols(0,MeshDimension-1).t()=\n"<<Coordinates.cols(0,MeshDimension-1).t();
        //cout<<"\ndN_by_dEps[ElementType][GaussPntr];\n"<<dN_by_dEps[ElementType][GaussPntr];
    }

     double Get_dS(int ElementType, int ElementNumber, int GaussPntr)
    {
        umat NodesAtElmntNmbr=Msh->ElmntPhysclGrpNodes[ElementType][PhysclGrpNum].row(ElementNumber);
        mat Coordinates=Msh->NodalCoordinates.rows(NodesAtElmntNmbr);
        mat F_rectangle=Coordinates.cols(0,originalDimension-1).t()*dN_by_dEps[ElementType][GaussPntr];
        //cout<<"dS =\n"<<F_rectangle.col(0)<<F_rectangle.col(1);
        return norm(cross(F_rectangle.col(0),F_rectangle.col(1)));
    }


private:
    int &originalDimension, PhysclGrpNum;
};

#endif // TRIALFUNCTIONNEUMANNSURFACE_HPP
