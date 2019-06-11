#ifndef TESTFUNCTIONGALERKIN_HPP
#define TESTFUNCTIONGALERKIN_HPP
#include "TrialFunction.hpp"
template<class GenericTrialFunction>
class TestFunctionGalerkin
{
public:
    int vectorLvl;
    TestFunctionGalerkin()
    {
        //---Do-Nothing
    }
    /// Must Initialize using this.
    TestFunctionGalerkin(GenericTrialFunction& u1)
    {
        u=&u1;
        vectorLvl=u1.vectorLvl;
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
