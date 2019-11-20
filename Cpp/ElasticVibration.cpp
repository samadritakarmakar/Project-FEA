#include "ProjectFEA.hpp"
using namespace arma;

class M_local: public LocalIntegrator<TrialFunction>
{
public:
    M_local(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, double& density):
        LocalIntegrator (a,u,v), rho(density)
    {
    }
    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v,
                     int thread)
    {
        return rho*a[thread].v(v)*a[thread].u(u)*a[thread].dX(u);
    }
    double & rho;
};


/// A Linear Elastic Model is defined here. This weak form is integrated over each element.
/// The Virtual 'weak_form' function defined in 'LocalIntegrator' is overloaded during runtime.
class K_local: public LocalIntegrator<TrialFunction>
{
public:
    K_local(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, sp_mat& C):
        LocalIntegrator (a,u,v), C(C)
    {
    }
    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v,
                     int thread)
    {
        return a[thread].inner(a[thread].sym_grad(v),C*a[thread].sym_grad(u))*a[thread].dX(u);
    }
    sp_mat& C;
};

/// A weak form for Body forces is defined here.
/// The Virtual 'weak_form_vector' function defined in 'LocalIntegrator' is overloaded during runtime.
class BodyForce: public LocalIntegrator<TrialFunction>
{
public:
    BodyForce(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, vec& b):
        LocalIntegrator (a,u,v), b(b)
    {
    }
    mat weak_form_vector(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        return a[thread].dot(v,b)*a[thread].dX(u);
    }
    vec& b;
};

/// A weak form for Traction forces over a surface is defined here.
/// The Virtual 'weak_form_vector' function defined in 'LocalIntegrator' is overloaded during runtime.
class TractionForceSurface: public LocalIntegrator <TrialFunctionNeumannSurface>
{
public:
    TractionForceSurface(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                            TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, vec& ForceVector):
        LocalIntegrator (a,u,v), ForceVector(ForceVector)
    {
    }

    mat weak_form_vector(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                         TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, int thread)
    {
        //double Fz=-4.0e3/(200*60);
        //double Fy=0.0;
        //double Fz=0.0;
        //double Fy=-0.046189;
        vec vctr=ForceVector;
        return a[thread].dot(v,vctr)*a[thread].dS(u);
    }
    vec &ForceVector;
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

void Build_Matrix_M(sp_mat &M, M_local& M_local_elm, SystemAssembler<M_local, TrialFunction>& M_build)
{
    M_build.SetMatrixSize(M);
    M_build.RunSystemAssembly(M_local_elm, M);
}

void Build_Matrix_K(sp_mat &K, K_local& K_local_elm, SystemAssembler<K_local, TrialFunction>& K_build)
{
    K_build.SetMatrixSize(K);
    K_build.RunSystemAssembly(K_local_elm, K);
}

class ForceVectorCalculator
{
public:
    ForceVectorCalculator(vec& MaxForce, double StartTime, double StopForceTime):
        StartTime(StartTime), StopForceTime(StopForceTime), MaxForce(MaxForce)
    {}
    double RampFunction(double currentTime)
    {
        return (currentTime-StartTime)/(StopForceTime-StartTime);
    }
    void Calculate_Force(vec& CurrentForce, double currentTime)
    {
        double Ramp;
        CurrentForce.set_size(MaxForce.n_rows);
        if(StopForceTime >= currentTime)
        {
            Ramp=RampFunction(currentTime);
            CurrentForce=Ramp*MaxForce;
        }
        else
        {
            CurrentForce.zeros();
        }
    }
    double StartTime, StopForceTime;
    vec& MaxForce;
};

void Build_Final_f(mat &f, mat& f1, mat& f2, vec& CurrentForceVector, double time, ForceVectorCalculator& CalcForce,
                   TractionForceSurface & f2_local, SystemAssembler<TractionForceSurface, TrialFunctionNeumannSurface>& f2_build)
{

    //Since current force vector is referenced when initalizing f2_local, the latest value of CurrentForceVector
    //shall be taken in consideration everytime RunSystemAssemblyVector is being executed.
    CalcForce.Calculate_Force(CurrentForceVector, time);
    //f2_build.SetVectorSize(f2);
    f2.zeros();
    f2_build.RunSystemAssemblyVector(f2_local, f2);
    f=f1+f2;
}

//void Build_Matrix_f(mat& f, )

int main(int argc, char *argv[])
{
    if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./ElasticVibration <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    libGmshReader::MeshReader Mesh(FileName, Dimension);

    double E=206.84e3, nu=0.3, density=7700.0* 1e-9, g=9.81; //density in kg/mm3.
    double beta=0.001, gamma=0.002; //the contributions to the damping matrix from M and K matrix respectively.
    int vectorLevel=3;
    int NeumannSurfaceGroup=1, DirichletSurfaceGroup=0, VolGroup=0;
    double Fx=0.0, Fz=0.0, Fy=-0.046189;
    vec MaxForceVector={Fx, Fy, Fz};
    vec CurrentForceVector;

    double StartTime=0.0, StopForceTime=1.0, delta_t=0.1, EndTime=5.0;
    std::vector<double> time(2), theta(2);
    time[0]=StartTime;
    time[1]=StartTime;

    theta[0]=0.5;
    theta[1]=0.5;

    sp_mat C;
    Set_C_3D(C, E, nu);
    std::vector<sp_mat> M_C_K(3);// K should be at index 0, C at 1 and M at 2.
    std::vector<mat> f(3);
    mat f1;

    std::vector<mat> Sol_of_u(2);

    TrialFunction u_ddot(Mesh, vectorLevel);
    TestFunctionGalerkin<TrialFunction> v_ddot(u_ddot);
    FormMultiThread<TrialFunction> a_ddot;
    M_local M_local_elm(a_ddot, u_ddot, v_ddot, density);
    SystemAssembler<M_local, TrialFunction> M_build(a_ddot, u_ddot, v_ddot);
    Build_Matrix_M(M_C_K[2], M_local_elm, M_build); //M Matrix of M ddot_x + C dot_x + K x =  f

    TrialFunction u(Mesh,vectorLevel);
    TestFunctionGalerkin<TrialFunction> v(u);
    FormMultiThread<TrialFunction> a;
    K_local K_local_elm(a, u, v, C);
    SystemAssembler<K_local, TrialFunction> K_build(a, u, v);
    Build_Matrix_K(M_C_K[0], K_local_elm, K_build); //K Matrix of M ddot_x + C dot_x + K x =  f
    K_build.SetMatrixSize(M_C_K[1]);
    M_C_K[1]=beta*M_C_K[0]+gamma*M_C_K[1]; //C Matrix of M ddot_x + C dot_x + K x =  f

    vec bodyForce={0,0,-1};
    bodyForce=density*g*bodyForce;
    BodyForce f1_local(a, u, v, bodyForce);
    SystemAssembler<BodyForce, TrialFunction> f1_build(a, u, v);
    f1_build.SetVectorSize(f1);
    f1_build.RunSystemAssemblyVector(f1_local, f1);

    ForceVectorCalculator CalcForce(MaxForceVector, StartTime, StopForceTime);

    TrialFunctionNeumannSurface u1(u, NeumannSurfaceGroup);
    TestFunctionGalerkin<TrialFunctionNeumannSurface> v1(u1);
    FormMultiThread<TrialFunctionNeumannSurface> a1;
    //Since current force vector is referenced when initalizing f2_local, the latest value of CurrentForceVector
    //shall be taken in consideration everytime RunSystemAssemblyVector is being executed.
    TractionForceSurface f2_local(a1, u1, v1, CurrentForceVector);
    SystemAssembler<TractionForceSurface, TrialFunctionNeumannSurface> f2_build(a1, u1, v1);
    mat f2;
    f1_build.SetVectorSize(f2);
    f1_build.SetVectorSize(f[0]);
    f1_build.SetVectorSize(f[1]);
    f1_build.SetVectorSize(f[2]);
    Build_Final_f(f[0], f1, f2, CurrentForceVector, time[1], CalcForce, f2_local, f2_build);
    time[1]+=delta_t;
    Build_Final_f(f[1], f1, f2, CurrentForceVector, time[1], CalcForce, f2_local, f2_build);
    time[1]+=delta_t;
    Build_Final_f(f[2], f1, f2, CurrentForceVector, time[1], CalcForce, f2_local, f2_build);

    umat boolDiricletNodes={1,1,1};
    DirichletBC DrchltBC(u1,0, boolDiricletNodes);
    DeclaredExprssn ExprssnZero(vectorLevel);
    DrchltBC.SetDirichletBCExpression(ExprssnZero);

    InitialBC InitBC(u, VolGroup, boolDiricletNodes);
    InitBC.SetInitialBCExpression(ExprssnZero);
    f1_build.SetVectorSize(Sol_of_u[0]);
    InitBC.ApplyInitialBC(Sol_of_u[0]);
    f1_build.SetVectorSize(Sol_of_u[1]);
    InitBC.ApplyInitialBC(Sol_of_u[1]);

    SingleStepAlgorithm S22_Instance(u, v, delta_t);

    GmshWriter WriteToGmsh(u, "ElstcVbrtn.msh");

    TrialFunction Volume(Mesh, 1);
    TestFunctionGalerkin<TrialFunction> v_Vol(Volume);
    FormMultiThread<TrialFunction> a_vol;
    FindSize Vol_Size(a_vol, Volume, v_Vol);
    SystemAssembler<FindSize, TrialFunction> VolAssmbly(a_vol, Volume, v_Vol);
    mat Vol;
    VolAssmbly.SetScalarSize(Vol);
    VolAssmbly.RunScalarIntegration(Vol_Size, Vol);


    int vectorLevelStrss=6;
    TrialFunction stress(Mesh,vectorLevelStrss);
    FormMultiThread<TrialFunction> a_stress;
    TestFunctionGalerkin<TrialFunction> v_Stress(stress);

    std::vector<umat> NodePositions_u;
    K_build.SetNodePositions_u(NodePositions_u);

    mat Strss;
    Stress strss(a_stress, u, v_Stress, C, Sol_of_u[0], Vol, NodePositions_u);
    SystemAssembler<Stress, TrialFunction> StrssAssmbly(a_stress, stress, v);
    //StrssAssmbly.SetVectorSize(Strss);
    StrssAssmbly.SetScalarSize(Strss);

    GmshWriter WriteStrss(stress, "ElstcVbrnStrss.msh");
    WriteStrss.viewName="Stress";
    WriteStrss.SetDataType_to_ElementData();

    for (int step=0; step<=EndTime/delta_t; step++)
    {
        cout<<"Currently at Step "<<step<<"\n";
        WriteToGmsh.WriteToGmsh(Sol_of_u[0], step, time[0]);
        WriteStrss.WriteToGmshSymmetric_3x3(Strss, step, time[0]);
        sp_mat A; mat b;
        S22_Instance.SetSizeOfMatrix(A);
        S22_Instance.SetSizeOfVector(b);
        S22_Instance.SingleStep_22(A, b, Sol_of_u, theta, M_C_K, f);
        DrchltBC.ApplyBC_Dynamic(Sol_of_u, A, b);
        mat alpha;
        spsolve(alpha, A, b);
        mat f_new;
        time[0]+=delta_t;
        time[1]+=delta_t;
        S22_Instance.update_Sol_of_u(Sol_of_u, alpha);
        f1_build.SetVectorSize(f_new);
        Build_Final_f(f_new, f1, f2, CurrentForceVector, time[1], CalcForce, f2_local, f2_build);
        S22_Instance.update_f(f, f_new);

        Strss.zeros();
        StrssAssmbly.RunScalarIntegration(strss, Strss);
    }
}
