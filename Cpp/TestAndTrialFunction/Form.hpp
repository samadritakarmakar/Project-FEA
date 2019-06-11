#ifndef FORM_HPP
#define FORM_HPP
#include "TrialFunction.hpp"
#include "TrialFunctionNeumannSurface.hpp"
#include "TrialFunctionNeumannLine.hpp"
#include "TestFunctionGalerkin.hpp"
#include <armadillo>

template<class GenericTrialFunction>
class Form//: public TrialFunction
{
public:
    sp_mat ResultingMat;
    mat ResultingVector;
    int ElementType;
    int ElementNumber;
    int GaussPntr;
    Form()
    {
        ElementType=0;
        ElementNumber=0;
        GaussPntr=0;
       /* u_down=std::vector<TrialFunction> (u.MeshDimension-1);
        v_down=std::vector<TestFunctionGalerkin> (u.MeshDimension-1);
        libGmshReader::MeshReader Mesh(u.Msh->ElementData::fileName,u.MeshDimension-1);
        u_down[0]=TrialFunction(Mesh, u.vectorLvl-1);
        v_down[0]=TestFunctionGalerkin(u_down[0]);

        cout<<mat(inner(grad(v_down[0]),grad(u_down[0]))); */
    }


    Form(GenericTrialFunction& u): u_Internal(&u)
    {
        ElementType=0;
        ElementNumber=0;
        GaussPntr=0;

       /* u_down=std::vector<TrialFunction> (u.MeshDimension-1);
        v_down=std::vector<TestFunctionGalerkin> (u.MeshDimension-1);
        libGmshReader::MeshReader Mesh(u.Msh->ElementData::fileName,u.MeshDimension-1);
        u_down[0]=TrialFunction(Mesh, u.vectorLvl-1);
        v_down[0]=TestFunctionGalerkin(u_down[0]);

        cout<<mat(inner(grad(v_down[0]),grad(u_down[0]))); */
    }


    void set_u_Internal(GenericTrialFunction& u)
    {
        u_Internal=&u;
    }

    void set_v_Internal(TestFunctionGalerkin<GenericTrialFunction>& v)
    {
        v_Internal=&v;
    }


     sp_mat u(GenericTrialFunction& u)
    {
        //set_u_Internal(u);
        return u.Get_u(ElementType,GaussPntr);
    }

     sp_mat v(TestFunctionGalerkin<GenericTrialFunction>& v)
    {
        //set_v_Internal(v);
        return v.Get_v(ElementType, GaussPntr);
    }

     mat x(GenericTrialFunction& u)
     {
         return u.Get_x(ElementType, ElementNumber, GaussPntr);
     }


     sp_mat grad(GenericTrialFunction& u)
    {
        //set_u_Internal(u);
        return u.Get_grad_u(ElementType, ElementNumber, GaussPntr);
    }

     sp_mat grad(TestFunctionGalerkin<GenericTrialFunction>& v)
    {
        //set_v_Internal(v);
        return v.Get_grad_v(ElementType, ElementNumber, GaussPntr);
    }

     sp_mat transpose_grad(GenericTrialFunction& u)
     {
         return u.GetTranspose_grad_u(ElementType, ElementNumber, GaussPntr);
     }

     sp_mat transpose_grad(TestFunctionGalerkin<GenericTrialFunction>& v)
     {
         return v.GetTranspose_grad_v(ElementType, ElementNumber, GaussPntr);
     }
     sp_mat curl(GenericTrialFunction& u)
    {
        //set_u_Internal(u);
        return u.Get_curl_u(ElementType, ElementNumber, GaussPntr);
    }

     sp_mat curl(TestFunctionGalerkin<GenericTrialFunction>& v)
    {
        return v.Get_curl_v(ElementType, ElementNumber, GaussPntr);
    }

     mat trace(TrialFunction& u)
     {
         return u.Get_trace_grad_u(ElementType, ElementNumber, GaussPntr);
     }

     sp_mat inner(TestFunctionGalerkin<GenericTrialFunction>& v, GenericTrialFunction& u)
    {
        //set_u_Internal(u);
        return v.Get_v(ElementType,GaussPntr)*u.Get_u(ElementType,GaussPntr);
    }

     sp_mat inner(TestFunctionGalerkin<GenericTrialFunction>& v, sp_mat u)
    {
        //set_v_Internal(v);
        return v.Get_v(ElementType,GaussPntr)*u;
    }

     sp_mat inner(sp_mat grad_v, sp_mat grad_u)
    {
        return grad_v*grad_u;
    }

     mat dot(TestFunctionGalerkin<GenericTrialFunction>& v, vec& vctr)
     {
         return v.Get_v(ElementType,GaussPntr)*vctr;
     }

     sp_mat dot(TestFunctionGalerkin<GenericTrialFunction>& v, sp_mat grad_u)
    {
        //set_v_Internal(v);
        return v.Get_v(ElementType,GaussPntr)*grad_u;
    }

     sp_mat dot(sp_mat grad_v, sp_mat grad_u)
    {
        return grad_v*grad_u;
    }

     sp_mat dot(vec a, sp_mat grad_u)
    {
        //return dot_vectrLvl_grad_u(a,grad_u);
        int& vectorLvl = u_Internal->vectorLvl;
        if (vectorLvl==1)
        {
            //mat aMatrx=repmat(a.t(),grad_u.n_rows,1);
            //return sp_mat(aMatrx%grad_u);
            mat aMatrx=a.rows(1,u_Internal->MeshDimension).t();
            //cout<<grad_u;
            return sp_mat(aMatrx*grad_u);
        }
        else
        {
            mat vctr=a(span(0,u_Internal->vectorLvl-1)).t();
            sp_mat vecMatrx(vctr.n_cols,vectorLvl*vctr.n_cols);
            for (int i=0; i<vectorLvl; i++)
            {
                //cout<<"col ="<<i*vectorLvl<<":"<<i*vectorLvl+vectorLvl-1;
                vecMatrx.cols(i*vectorLvl,i*vectorLvl+vectorLvl-1)=sp_mat(diagmat(vctr.t()));
            }
            return vecMatrx*grad_u;
        }

    }


     double dX(GenericTrialFunction u)
    {
        mat F;
        //set_u_Internal(u);
        //cout<<"Gauss pt is at "<<GaussPntr<<"\n";
        u.Get_F(ElementType, ElementNumber, GaussPntr, F);
        return u_Internal->WeightAt(ElementType, GaussPntr)*norm(F);
    }


     double dL(GenericTrialFunction u)
    {
        //set_u_Internal(u);
        // cout<<"Gauss pt is at "<<GaussPntr<<"\n";
        return u_Internal->WeightAt(ElementType, GaussPntr)*u.Get_dL(ElementType, ElementNumber, GaussPntr);
    }


     double dS(GenericTrialFunction u)
    {
        //set_u_Internal(u);
         //cout<<"Gauss pt is at "<<GaussPntr<<"\n";
        return u_Internal->WeightAt(ElementType, GaussPntr)* u.Get_dS(ElementType, ElementNumber, GaussPntr);
    }

     sp_mat operator=(sp_mat Expression)
    {
        ResultingMat=Expression;
        return ResultingMat;
    }

     mat operator=(mat Expression)
     {
         ResultingVector=Expression;
         return ResultingVector;
     }

     void NextElementType()
    {
        ElementType++;
    }

     void NextElementNumber()
    {
        ElementNumber++;
    }

     void NextGaussPntr()
    {
        GaussPntr++;
    }

protected:
    GenericTrialFunction *u_Internal;
    TestFunctionGalerkin<GenericTrialFunction> *v_Internal;
};



#endif // FORM_HPP
