#include <iostream>
#include "ProjectFEA.hpp"

using namespace arma;
class new_LocalIntegrator: public LocalIntegrator<TrialFunction>
{
public:
    new_LocalIntegrator(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator (a,u,v)
    {
    }
    sp_mat weak_form(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v)
    {
        //return a.inner(a.curl(v),a.curl(u))*a.dX(u);
        sp_mat C;
        double E=200.0e3;
        Set_C_3D(C,E);
        return a.inner(a.sym_grad(v),C*a.sym_grad(u))*a.dX(u);
    }
};

class new_LocalIntegrator2: public LocalIntegrator<TrialFunction>
{
public:
    new_LocalIntegrator2(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator (a,u,v)
    {
    }
    mat weak_form_vector(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v)
    {
        vec b;
        b<<0<<endr<<0<<endr<<0<<endr;
        return a.dot(v,b)*a.dX(u);
    }
};

class new_Neu_Surf_LclIntgrtr: public LocalIntegrator <TrialFunctionNeumannSurface>
{
public:
    new_Neu_Surf_LclIntgrtr(Form<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                            TestFunctionGalerkin<TrialFunctionNeumannSurface>& v):
        LocalIntegrator (a,u,v)
    {
    }

    mat weak_form_vector(Form<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                         TestFunctionGalerkin<TrialFunctionNeumannSurface>& v)
    {
        double Fz=-4.0e3/(200*60);
        vec vctr={0,0,Fz};
        return a.dot(v,vctr)*a.dS(u);
    }

};


class new_Neu_Line_LclIntgrtr: public LocalIntegrator<TrialFunctionNeumannLine>
{
public:
    new_Neu_Line_LclIntgrtr(Form<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
                            TestFunctionGalerkin<TrialFunctionNeumannLine>& v):
        LocalIntegrator (a,u,v)
    {
    }
    mat weak_form_vector(Form<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
                     TestFunctionGalerkin<TrialFunctionNeumannLine>& v)
    {
        //vec vctr=vec(a.x(u));
        vec vctr={0,0,-4.0e3/2};
        return a.dot(v,vctr)*a.dL(u);
    }
};

class DeclaredExprssn : public Expression
{
public:
    DeclaredExprssn (int vectorLevel): Expression (vectorLevel)
    {
    }

    vec Eval(vec& x)
    {
        vec value={0,0,0};
        return value;
    }
};

int main(int argc, char *argv[])
{
    if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./ProjectTestFunction <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    //cout<<"File Name= "<<FileName<<"\n";
    //cout<<"Dimension= "<<Dimension<<"\n";
    libGmshReader::MeshReader Mesh(FileName, Dimension);

    int vectorLevel=3;
    TrialFunction u(Mesh,vectorLevel);
    TestFunctionGalerkin<TrialFunction> v(u);
    Form<TrialFunction> a;
    new_LocalIntegrator lcl_intgrt(a,u,v);
    SystemAssembler<new_LocalIntegrator, TrialFunction> systmAssmbly(a,u,v);
    //systmAssmbly.SetLocalIntegrator(lcl_intgrt);
    //sp_mat A;
    VariableMatrix A(1);
    systmAssmbly.SetMatrixSize(A.Matrix[0][0]);
    systmAssmbly.RunSystemAssembly(lcl_intgrt, A.Matrix[0][0]);
    //cout<<"zeros rows="<<all(A.Matrix[0][0]*ones(A.Matrix[0][0].n_cols,1)==0,1);
    //cout<<"zeros cols="<<all(A.Matrix[0][0].t()*ones(A.Matrix[0][0].n_cols,1)==0,1);

    Form<TrialFunction> a2;
    new_LocalIntegrator2 lcl_intgrt2(a2,u,v);
    SystemAssembler<new_LocalIntegrator2, TrialFunction> systmAssmbly2(a2,u,v);
    //mat b;
    VariableVector b(1);
    systmAssmbly2.SetVectorSize(b.Vector[0]);
    systmAssmbly2.RunSystemAssemblyVector(lcl_intgrt2,b.Vector[0]);

    Form<TrialFunctionNeumannSurface> a3;
    TrialFunctionNeumannSurface u_surf(u,1);
    //TestFunctionGalerkin<TrialFunctionNeumannSurface> v_surf(u_surf);
    //new_Neu_Surf_LclIntgrtr lclintgtr3(a3,u_surf, v_surf);
    //SystemAssembler<new_Neu_Surf_LclIntgrtr, TrialFunctionNeumannSurface> systmAssmbly3(a3, u_surf, v_surf);
    //systmAssmbly3.RunSystemAssemblyVector(lclintgtr3, b.Vector[0]);

    Form<TrialFunctionNeumannLine> a4;
    TrialFunctionNeumannLine u_line(u,0);
    TestFunctionGalerkin<TrialFunctionNeumannLine> v_line(u_line);
    new_Neu_Line_LclIntgrtr lcl_intgrt4(a4,u_line,v_line);
    SystemAssembler<new_Neu_Line_LclIntgrtr, TrialFunctionNeumannLine> systmAssmbly4(a4,u_line, v_line);
    systmAssmbly4.RunSystemAssemblyVector(lcl_intgrt4,b.Vector[0]);

    //cout<<b.Vector[0];

    umat boolDiricletNodes={1,1,1};
    DirichletBC DrchltBC(u_surf,0, boolDiricletNodes);
    DeclaredExprssn Dirich(vectorLevel);
    DrchltBC.SetDirichletBC(Dirich);
    DrchltBC.ApplyBC(A.Matrix[0][0],b.Vector[0]);
    //cout<<b.Vector[0];
    mat X=spsolve(A.Matrix[0][0],b.Vector[0]);
    //cout<<X.t();
    GmshWriter Write(u, "ElasticLine.pos");
    Write.WriteToGmsh(X);

    cout<<"Done!!!\n";
    return 0;
}
