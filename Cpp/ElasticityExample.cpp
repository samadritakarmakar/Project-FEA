#include <iostream>
#include "ProjectFEA.hpp"

using namespace arma;
/// A new LinearElastic Model is defined here. This weak form is integrated over each element.
/// The Virtual 'weak_form' function defined in 'LocalIntegrator' is overloaded during runtime.
class LinearElastic: public LocalIntegrator<TrialFunction>
{
public:
    LinearElastic(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, double E, double nu):
        LocalIntegrator (a,u,v), E(E), nu(nu)
    {
    }
    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v,
                     int thread)
    {
        sp_mat C;
        Set_C_3D(C,E,nu);
        return a[thread].inner(a[thread].sym_grad(v),C*a[thread].sym_grad(u))*a[thread].dX(u);
    }
    double E,nu;
};

/// A weak form for Body forces is defined here.
/// The Virtual 'weak_form_vector' function defined in 'LocalIntegrator' is overloaded during runtime.
class BodyForce: public LocalIntegrator<TrialFunction>
{
public:
    BodyForce(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator (a,u,v)
    {
    }
    mat weak_form_vector(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        vec b;
        b<<0<<endr<<0<<endr<<0<<endr;
        return a[thread].dot(v,b)*a[thread].dX(u);
    }
};

/// A weak form for Traction forces over a surface is defined here.
/// The Virtual 'weak_form_vector' function defined in 'LocalIntegrator' is overloaded during runtime.
class TractionForceSurface: public LocalIntegrator <TrialFunctionNeumannSurface>
{
public:
    TractionForceSurface(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                            TestFunctionGalerkin<TrialFunctionNeumannSurface>& v):
        LocalIntegrator (a,u,v)
    {
    }

    mat weak_form_vector(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                         TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, int thread)
    {
        //double Fz=-4.0e3/(200*60);
        //double Fy=0.0;
        double Fz=0.0;
        double Fy=-0.046189;
        vec vctr={0, Fy, Fz};
        return a[thread].dot(v,vctr)*a[thread].dS(u);
    }

};

/// A weak form for Traction forces over a Line is defined here.
/// The Virtual 'weak_form_vector' function defined in 'LocalIntegrator' is overloaded during runtime.
class TractionForceLine: public LocalIntegrator<TrialFunctionNeumannLine>
{
public:
    TractionForceLine(Form<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
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

class Stress : public LocalIntegrator<TrialFunction>
{
public:
    Stress(FormMultiThread<TrialFunction>& a, TrialFunction& u,
           TestFunctionGalerkin<TrialFunction>& v, sp_mat&C, mat& X, mat& vol, std::vector<umat>& NodePositions_u):
        LocalIntegrator (a, u, v), Disp(X), vol(vol), C(C), NodePositions_u(NodePositions_u)
    {}

    mat scalar_integration(FormMultiThread<TrialFunction>& a, TrialFunction& u, int thread)
    {
        umat positions=NodePositions_u[a[thread].ElementType].row(a[thread].ElementNumber);
        mat u_local=Disp.rows(positions);
        double element_volume=vol(a[thread].ElementNumber);
        return C*a[thread].sym_grad(u)*u_local*a[thread].dX(u)/element_volume;
    }
    mat& Disp, & vol;
    sp_mat& C;
    std::vector<umat>& NodePositions_u;
};

class FindSize : public LocalIntegrator<TrialFunction>
{
public:
    FindSize(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator<TrialFunction> (a, u, v){}
    mat scalar_integration(FormMultiThread<TrialFunction>& a, TrialFunction& u, int thread)
    {
        mat dx;
        dx<<a[thread].dX(u)<<endr;
        return dx;
    }
};

/*class Stress : public LocalIntegrator<TrialFunction>
{
public:
    Stress(FormMultiThread<TrialFunction>& a, TrialFunction& u,
           TestFunctionGalerkin<TrialFunction>& v): LocalIntegrator (a, u, v){}

    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u,
                     TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        return a[thread].v(v)*a[thread].u(u)*a[thread].dX(u);
    }
};

class StressMat : public LocalIntegrator<TrialFunction>
{
public:
    StressMat(FormMultiThread<TrialFunction>& a, TrialFunction& u,
           TestFunctionGalerkin<TrialFunction>& v, sp_mat& C):
        LocalIntegrator (a,u,v), C_Internal(C)
    {
        inner_C_grad_u=std::vector<sp_mat>(a.GetNumOfThreads());
    }

    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u,
                     TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        inner_C_grad_u[thread]=(C_Internal*a[thread].sym_grad(u));
        return a[thread].dot(v,inner_C_grad_u[thread])*a[thread].dX(u);
    }
private:
    sp_mat C_Internal;
    std::vector<sp_mat> inner_C_grad_u;
};*/
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
        std::cout<<"Usage: ./ElasticityEx <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    libGmshReader::MeshReader Mesh(FileName, Dimension);

    /// Vector Level Defines the Degree of Freedom per Node.
    int vectorLevel=3;

    /// Declaration of Quantity to be Calculated.
    /// Here it is Displacement.
    TrialFunction u(Mesh,vectorLevel);
    /// Here the Galerkin Test Fucntion is set.
    TestFunctionGalerkin<TrialFunction> v(u);
    /// Form is an interface Trial Function Declarations
    /// and the SystemAssembly.
    //Form<TrialFunction> a;
    FormMultiThread<TrialFunction> a;
    double E=206.84e3, nu=0.3;
    //E=200e3;
    /// Here an instance of the local intergrator is being declared.
    /// This is used to integrate over each element
    LinearElastic lcl_intgrt(a,u,v, E, nu);
    /// An instace of system Assembler is declared.
    SystemAssembler<LinearElastic, TrialFunction> systmAssmbly(a,u,v);
    //systmAssmbly.SetLocalIntegrator(lcl_intgrt);
    //sp_mat A;
    /// A matrix is defined. Alternatively, a sparse matrix using 'sp_mat'
    /// can also be defined.
    VariableMatrix A(1);
    /// The size matrix A is declared.
    systmAssmbly.SetMatrixSize(A.Matrix[0][0]);
    /// Here all the elements are integrated over and filled into the Matrix 'A'.
    systmAssmbly.RunSystemAssembly(lcl_intgrt, A.Matrix[0][0]);


    FormMultiThread<TrialFunction> a2;
    BodyForce lcl_intgrt2(a2,u,v);
    SystemAssembler<BodyForce, TrialFunction> systmAssmbly2(a2,u,v);
    //mat b;
    VariableVector b(1);
    /// The size of vector b is set here.
    systmAssmbly2.SetVectorSize(b.Vector[0]);
    /// Here all the elements are integrated over and filled into the vector 'b'.
    systmAssmbly2.RunSystemAssemblyVector(lcl_intgrt2,b.Vector[0]);

    FormMultiThread<TrialFunctionNeumannSurface> a3;
    /// Here a Traction force over the surface is being declared.
    /// The first argument specifies the function 'u' that is to be solved for.
    /// The second argument specifies the Index of the Physical Group defined in Gmsh.
    /// The Physical Group Index here does not depend on the Physical Group number of Gmsh
    /// but is rather an Index of where the group was defined 1st 2nd or 3rd (in here, 0, 1,or 2)
    TrialFunctionNeumannSurface u_surf(u,1);
    /// Here the Galerkin Test Function is set.
    TestFunctionGalerkin<TrialFunctionNeumannSurface> v_surf(u_surf);
    /// An Instance For calculating the Traction forces on the surface is defined here.
    TractionForceSurface lclintgtr3(a3,u_surf, v_surf);
    SystemAssembler<TractionForceSurface, TrialFunctionNeumannSurface> systmAssmbly3(a3, u_surf, v_surf);
    /// Here all the elements are integrated over and filled into the vector 'b'.
    systmAssmbly3.RunSystemAssemblyVector(lclintgtr3, b.Vector[0]);

    //Form<TrialFunctionNeumannLine> a4;
    //TrialFunctionNeumannLine u_line(u,0);
    //TestFunctionGalerkin<TrialFunctionNeumannLine> v_line(u_line);
    //TractionForceLine lcl_intgrt4(a4,u_line,v_line);
    //SystemAssembler<TractionForceLine, TrialFunctionNeumannLine> systmAssmbly4(a4,u_line, v_line);
    //systmAssmbly4.RunSystemAssemblyVector(lcl_intgrt4,b.Vector[0]);

    //cout<<b.Vector[0];

    /// boolDiricletNodes specifies over which of the degrees of freedom of a node the Boundary Condition has to be applied .
    umat boolDiricletNodes={1,1,1};
    /// Here an Instance of DirichletBC is defined.
    DirichletBC DrchltBC(u_surf,0, boolDiricletNodes);
    /// The degree of freedom for each node for the DirichletBC is defined.
    DeclaredExprssn Dirich(vectorLevel);
    /// The expression evaluated to be used at each node to set the DirichletBC is set here.
    DrchltBC.SetDirichletBCExpression(Dirich);
    /// The dirchlet BC is set here.
    DrchltBC.ApplyBC(A.Matrix[0][0],b.Vector[0]);
    /// The problem is solved here.
    mat X;
    spsolve(X, A.Matrix[0][0],b.Vector[0]);

    GmshWriter Write(u, "ElasticEx.pos");
    Write.viewName="Displacement";
    Write.WriteToGmsh(X);
    //Write.WriteToGmsh(X, "Disp");

    TrialFunction Volume(Mesh, 1);
    TestFunctionGalerkin<TrialFunction> v_Vol(Volume);
    FormMultiThread<TrialFunction> a_vol;
    FindSize Vol_Size(a_vol, Volume, v_Vol);
    SystemAssembler<FindSize, TrialFunction> VolAssmbly(a_vol, Volume, v_Vol);
    mat Vol;
    VolAssmbly.SetScalarSize(Vol);
    VolAssmbly.RunScalarIntegration(Vol_Size, Vol);


    //Setting Mesh to order 1 for Post process analysis
    vectorLevel=6;
    TrialFunction stress(Mesh,vectorLevel);
    FormMultiThread<TrialFunction> a_stress;
    TestFunctionGalerkin<TrialFunction> v_Stress(stress);

    std::vector<umat> NodePositions_u;
    systmAssmbly.SetNodePositions_u(NodePositions_u);

    sp_mat C;
    Set_C_3D(C, E);

    mat Strss;
    Stress strss(a_stress, u, v_Stress, C, X, Vol, NodePositions_u);
    SystemAssembler<Stress, TrialFunction> StrssAssmbly(a_stress, stress, v);
    StrssAssmbly.SetVectorSize(Strss);
    StrssAssmbly.RunScalarIntegration(strss, Strss);

    GmshWriter WriteStrss(stress, "ElstcStrss.pos");
    WriteStrss.viewName="Stress";
    WriteStrss.SetDataType_to_ElementData();
    WriteStrss.WriteToGmshSymmetric_3x3(Strss);

    cout<<"Done!!!\n";
    return 0;
}
