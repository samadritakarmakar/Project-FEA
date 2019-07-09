#include <iostream>
#include "ProjectFEA.hpp"


/// A new Model is defined here. This weak form is integrated over each element.
/// The Virtual 'weak_form' function defined in 'LocalIntegrator' is overloaded during runtime.
class new_LocalIntegrator: public LocalIntegrator<TrialFunction>
{
public:
    new_LocalIntegrator(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator (a,u,v)
    {
    }
    sp_mat weak_form(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v)
    {
        return a.inner(a.grad(v),a.grad(u))*a.dX(u);
    }
};

/// The Virtual 'weak_form_vector' function defined in 'LocalIntegrator' is overloaded during runtime.
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

/// The Virtual 'weak_form_vector' function defined in 'LocalIntegrator' is overloaded during runtime.
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
        vec vctr={1,1,1};
        return a.dot(v,vctr)*a.dS(u);
    }

};

/// The Virtual 'weak_form_vector' function defined in 'LocalIntegrator' is overloaded during runtime.
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
        vec vctr={1,1,1};
        return a.dot(v,vctr)*a.dL(u);
    }
};

///This class over here through its overloaded virtual function declares the values of Dirichlet Nodes.
/// The virtual function 'Eval' is evaluated at each node to find the value of Dirichlet Condtion at that node.
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
        std::cout<<"Usage: ./Poisson <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    //cout<<"File Name= "<<FileName<<"\n";
    //cout<<"Dimension= "<<Dimension<<"\n";
    libGmshReader::MeshReader Mesh(FileName, Dimension);

    int vectorLevel=3;
    /// Declaration of Quantity to be Calculated.
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

    Form<TrialFunction> a2;
    new_LocalIntegrator2 lcl_intgrt2(a2,u,v);
    SystemAssembler<new_LocalIntegrator2, TrialFunction> systmAssmbly2(a2,u,v);
    //mat b;
    VariableVector b(1);
    systmAssmbly2.SetVectorSize(b.Vector[0]);
    systmAssmbly2.RunSystemAssemblyVector(lcl_intgrt2,b.Vector[0]);

    Form<TrialFunctionNeumannSurface> a3;
    TrialFunctionNeumannSurface u_surf(u,0);
    TestFunctionGalerkin<TrialFunctionNeumannSurface> v_surf(u_surf);
    new_Neu_Surf_LclIntgrtr lclintgtr3(a3,u_surf, v_surf);
    SystemAssembler<new_Neu_Surf_LclIntgrtr, TrialFunctionNeumannSurface> systmAssmbly3(a3, u_surf, v_surf);
    systmAssmbly3.RunSystemAssemblyVector(lclintgtr3, b.Vector[0]);

    //Form<TrialFunctionNeumannLine> a4;
    //TrialFunctionNeumannLine u_line(u,0);
    /*TestFunctionGalerkin<TrialFunctionNeumannLine> v_line(u_line);
    new_Neu_Line_LclIntgrtr lcl_intgrt4(a4,u_line,v_line);
    SystemAssembler<new_Neu_Line_LclIntgrtr, TrialFunctionNeumannLine> systmAssmbly4(a4,u_line, v_line);
    systmAssmbly4.RunSystemAssemblyVector(lcl_intgrt4,b.Vector[0]);*/

    //cout<<b.Vector[0];

    umat boolDiricletNodes={1,1,1};
    DirichletBC DrchltBC(u_surf,1, boolDiricletNodes);
    DeclaredExprssn Dirich(vectorLevel);
    DrchltBC.SetDirichletBCExpression(Dirich);
    DrchltBC.ApplyBC(A.Matrix[0][0],b.Vector[0]);
    //cout<<b.Vector[0];
    mat X=spsolve(A.Matrix[0][0],b.Vector[0]);
    //cout<<X;
    GmshWriter Write(u, "poisson.pos");
    Write.WriteToGmsh(X);
    //Write.WriteToGmsh(X,"Disp");

    cout<<"Done!!!\n";
    return 0;
}
