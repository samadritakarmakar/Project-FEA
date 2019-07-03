#include "ProjectFEA.hpp"
using namespace arma;

class LinearElastic: public LocalIntegrator<TrialFunction>
{
public:
    LinearElastic(Form<TrialFunction>& a, TrialFunction& u,
                  TestFunctionGalerkin<TrialFunction>& v, sp_mat& C):
        LocalIntegrator (a,u,v), C_Internal(C)
    {
    }
    sp_mat weak_form(Form<TrialFunction>& a, TrialFunction& u,
                     TestFunctionGalerkin<TrialFunction>& v)
    {
        //cout<<mat(a[thread].inner(a[thread].sym_grad(v),C_Internal*a[thread].sym_grad(u))*a[thread].dX(u))<<"\n";
        //cout<<mat(a.inner(a.sym_grad(v),C_Internal*a.sym_grad(u))*a.dX(u));
        return a.inner(a.sym_grad(v),C_Internal*a.sym_grad(u))*a.dX(u);
        //return a[thread].inner(a[thread].sym_grad(v),C_Internal*a[thread].sym_grad(u))*a[thread].dX(u);
    }

private:
    sp_mat& C_Internal;
};

class LinearElastic_BodyForce: public LocalIntegrator<TrialFunction>
{
public:
    LinearElastic_BodyForce(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, vec& BodyForce):
        BodyForce_Internal(BodyForce),
        LocalIntegrator (a,u,v)
    {
    }
    mat weak_form_vector(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v)
    {
        return a.dot(v,BodyForce_Internal)*a.dX(u);
    }
private:
    vec &BodyForce_Internal;
};

class LinearElastic_TractionArea: public LocalIntegrator <TrialFunctionNeumannSurface>
{
public:
    LinearElastic_TractionArea(Form<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                            TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, vec& TractionForce):
        TractionForce_Internal(TractionForce),
        LocalIntegrator (a,u,v)
    {
    }

    mat weak_form_vector(Form<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                         TestFunctionGalerkin<TrialFunctionNeumannSurface>& v)
    {
        return a.dot(v,TractionForce_Internal)*a.dS(u);
    }
private:
vec& TractionForce_Internal;
};


class LinearElastic_TractionLine: public LocalIntegrator<TrialFunctionNeumannLine>
{
public:
    LinearElastic_TractionLine(Form<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
                            TestFunctionGalerkin<TrialFunctionNeumannLine>& v,  vec& TractionForce):
        TractionForce_Internal(TractionForce),
        LocalIntegrator (a,u,v)
    {
    }
    mat weak_form_vector(Form<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
                     TestFunctionGalerkin<TrialFunctionNeumannLine>& v)
    {
        return a.dot(v,TractionForce_Internal)*a.dL(u);
    }
private:
vec& TractionForce_Internal;
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


int main()
{
    int vectorLevel=3;
    libGmshReader::MeshReader Mesh("OneElmntGeom/Hexahedral.msh",vectorLevel);
    TrialFunction u(Mesh, vectorLevel);
    TestFunctionGalerkin<TrialFunction> v(u);
    //Form<TrialFunction> a(u);
    Form<TrialFunction> a;
    VariableMatrix A(1);
    double E=200.0e3;
    sp_mat C;
    Set_C_3D(C, E);
    LinearElastic intgrt(a,u,v,C);
    SystemAssembler<LinearElastic, TrialFunction> systmAssmly(a,u,v);
    systmAssmly.SetMatrixSize(A.Matrix[0][0]);
    systmAssmly.RunSystemAssembly(intgrt,A.Matrix[0][0]);
    //cout<<mat(A.Matrix[0][0]);

    vec body_force;
    body_force={0,0,0};
    Form<TrialFunction> a2;
    LinearElastic_BodyForce lcl_intgrt2(a2,u,v, body_force);
    SystemAssembler<LinearElastic_BodyForce, TrialFunction> systmAssmbly2(a2,u,v);
    VariableVector b(1);
    systmAssmbly2.SetVectorSize(b.Vector[0]);
    systmAssmbly2.RunSystemAssemblyVector(lcl_intgrt2,b.Vector[0]);

    double Fx=4.0e3/(2*2);
    vec TractionForce={Fx,0,0};
    Form<TrialFunctionNeumannSurface> a3;
    TrialFunctionNeumannSurface u_surf(u,0);
    TestFunctionGalerkin<TrialFunctionNeumannSurface> v_surf(u_surf);
    LinearElastic_TractionArea lclintgtr3(a3,u_surf, v_surf, TractionForce);
    SystemAssembler<LinearElastic_TractionArea, TrialFunctionNeumannSurface> systmAssmbly3(a3, u_surf, v_surf);
    systmAssmbly3.RunSystemAssemblyVector(lclintgtr3, b.Vector[0]);

    //Form<TrialFunctionNeumannLine> a4;
    TrialFunctionNeumannLine u_line(u,0);
    //TestFunctionGalerkin<TrialFunctionNeumannLine> v_line(u_line);
    //LinearElastic_TractionLine lcl_intgrt4(a4,u_line,v_line);
    //SystemAssembler<LinearElastic_TractionLine, TrialFunctionNeumannLine> systmAssmbly4(a4,u_line, v_line);
    //systmAssmbly4.RunSystemAssemblyVector(lcl_intgrt4,b.Vector[0]);
    //cout<<b.Vector[0];

    umat boolDiricletNodes={1,1,1};
    DirichletBC DrchltBC(u_line,0, boolDiricletNodes);
    DeclaredExprssn Dirich(vectorLevel);
    DrchltBC.SetDirichletBCExpression(Dirich);
    DrchltBC.ApplyBC(A.Matrix[0][0],b.Vector[0]);
    //cout<<b.Vector[0];

    boolDiricletNodes={1,0,0};
    DirichletBC DrchltBC2(u_line,1, boolDiricletNodes);
    DrchltBC2.SetDirichletBCExpression(Dirich);
    DrchltBC2.ApplyBC(A.Matrix[0][0],b.Vector[0]);

    mat(A.Matrix[0][0]).save("A", arma_ascii);
    b.Vector[0].save("b", arma_ascii);
    mat X=spsolve(A.Matrix[0][0],b.Vector[0]);
    //cout<<X.t();
    GmshWriter Write(u, "ElmntVrfy.pos");
    Write.WriteToGmsh(X);
    return 0;
}
