#ifndef TESTFUNCTIONGALERKIN_HPP
#define TESTFUNCTIONGALERKIN_HPP
#include "TrialFunction.hpp"
template<class GenericTrialFunction>
class TestFunctionGalerkin
{
public:
    int MeshDimension;
    int vectorLvl;
    int originalVctrLvl;
    int numOfNodes;
    std::vector<int> numOfGaussPoints;
    std::vector<umat> ElmntNodes;
    TestFunctionGalerkin()
    {
        //---Do-Nothing
    }
    /// Must Initialize using this.
    TestFunctionGalerkin(GenericTrialFunction& _u)
    {
        u=&_u;
        MeshDimension=_u.MeshDimension;
        vectorLvl=_u.vectorLvl;
        originalVctrLvl=_u.originalVctrLvl;
        numOfNodes=_u.Msh->NodeTag.n_rows;
        numOfGaussPoints=std::vector<int>(_u.NoOfElementTypes);
        ElmntNodes=std::vector<umat>(_u.NoOfElementTypes);
        for (int ElmntTyp=0; ElmntTyp<_u.NoOfElementTypes; ElmntTyp++)
        {
            numOfGaussPoints[ElmntTyp]=_u.GetNumberOfGaussPoints(ElmntTyp);
            ElmntNodes[ElmntTyp]=_u.ElmntNodes[ElmntTyp];

        }
        //cout<<"v =\n"<<mat(Get_v(0,0))<<"\n";
        //cout<<"Grad_v\n"<<mat(Get_grad_v(0,0,0))<<"\n";
    }
    /// Return Shape Function of v in Matrix form
     sp_mat Get_v(int ElementType, int GaussPntr)
    {
        return u->Get_u(ElementType, GaussPntr).t();
    }
    /// Return Shape Function of grad of v in Matrix form
     sp_mat Get_grad_v(int ElementType, int ElementNumber, int GaussPntr)
    {
        return u->Get_grad_u(ElementType, ElementNumber, GaussPntr).t();
    }

     sp_mat Get_Sym_grad_v(int ElementType, int ElementNumber, int GaussPntr)
     {
         return u->Get_Sym_grad_u(ElementType, ElementNumber, GaussPntr).t();
     }

     sp_mat GetTranspose_grad_v(int ElementType,int ElementNumber,int GaussPntr)
     {
         return u->GetTranspose_grad_u(ElementType, ElementNumber, GaussPntr).t();
     }

     mat Get_trace_grad_v(int ElementType,int ElementNumber,int GaussPntr)
     {
         return u->Get_trace_grad_u(ElementType, ElementNumber, GaussPntr).t();
     }

     sp_mat Get_curl_v(int ElementType, int ElementNumber, int GaussPntr)
    {
        return u->Get_curl_u(ElementType, ElementNumber, GaussPntr).t();
    }

    void Get_F(int ElementType, int ElementNumber, int GaussPntr, mat& F)
    {
        u->Get_F(ElementType, ElementNumber, GaussPntr, F);
    }

private:
    GenericTrialFunction *u;
};

#endif // TESTFUNCTION_HPP
