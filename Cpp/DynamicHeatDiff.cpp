#include "ProjectFEA.hpp"
class Matrix_C: public LocalIntegrator<TrialFunction>
{
public:
    Matrix_C(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v,
             double& heatCapacity, double& density):
        LocalIntegrator (a,u,v), c(heatCapacity), rho(density)
    {
    }
    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        return rho*c*a[thread].v(v)*a[thread].u(u)*a[thread].dX(u);
    }
    double  &c; double & rho;
};

class Matrix_K1: public LocalIntegrator<TrialFunction>
{
public:
    Matrix_K1(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, double& thermal_conductivity):
        LocalIntegrator (a,u,v), k(thermal_conductivity)
    {
    }
    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        return k*a[thread].inner(a[thread].grad(v),a[thread].grad(u))*a[thread].dX(u);
    }
    double & k;
};

class Matrix_K2: public LocalIntegrator<TrialFunctionNeumannSurface>
{
public:
    Matrix_K2(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
              TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, double& ht_cnvctn_coeff):
        LocalIntegrator (a,u,v), h(ht_cnvctn_coeff)
    {
    }
    sp_mat weak_form(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                     TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, int thread)
    {
        return h*a[thread].v(v)*a[thread].u(u)*a[thread].dS(u);
    }
    double & h;
};




class Matrix_f_Source: public LocalIntegrator<TrialFunction>
{
public:
    Matrix_f_Source(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, double &source):
        LocalIntegrator (a,u,v), s(source)
    {
    }
    mat weak_form_vector(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v)
    {
        return mat(a.v(v)*s*a.dX(u));
    }
    double & s;
};

class Matrix_f_Neumann: public LocalIntegrator<TrialFunctionNeumannSurface>
{
public:
    Matrix_f_Neumann(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                     TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, mat &heat_normal,
                     double & conv_coeff, mat &u_ambient, std::vector<umat>& Positions):
        LocalIntegrator (a,u,v), u_a(u_ambient), h_conv(conv_coeff), q_n(heat_normal),
        MatrixPositions(Positions)
    {
    }
    mat weak_form_vector(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                         TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, int thread)
    {
        umat Position=MatrixPositions[a[thread].ElementType].row(a[thread].ElementNumber);
        return (a[thread].v(v)*a[thread].u(u)*(-q_n(Position)+h_conv*u_a(Position))*a[thread].dS(u));
    }
    mat & u_a;
    double &h_conv;
    mat & q_n;
    std::vector<umat>& MatrixPositions;
};

class InitialBC_Expression1 : public Expression
{
public:
    InitialBC_Expression1(double& temperature, double& time, int vectorLevel):
        Expression(vectorLevel), temp(temperature), time(time)
    {
    }

    vec Eval(vec& x)
    {
        vec t={temp};
        return t;
    }
    double &temp, &time;
};

class InitialBC_Expression2 : public Expression
{
public:
    InitialBC_Expression2(double& normal_heat_loss_surf, double & time, int vectorLevel):
        Expression(vectorLevel), q_n_bar(normal_heat_loss_surf), time(time)
    {
    }

    vec Eval(vec& x)
    {
        vec q_n={q_n_bar};
        return q_n;
    }
    double &q_n_bar, &time;
};

class DirichletBC_Expression : public Expression
{
public:
    DirichletBC_Expression(double AppliedTemperature, int vectorLevel): Expression(vectorLevel), AppliedTemperature(AppliedTemperature)
    {}

    vec Eval(vec& x)
    {
        vec AppliedTemp={AppliedTemperature};
        return AppliedTemp;
    }
    double AppliedTemperature;
};

void BuildMatrix_C(sp_mat& C, Matrix_C& C_local, SystemAssembler<Matrix_C, TrialFunction>& C_build)
{
    C_build.SetMatrixSize(C);
    C_build.RunSystemAssembly(C_local, C);
    cout<<"Matrix C built!\n";
}

void BuildMatrix_K(sp_mat& K, Matrix_K1& K1_local, SystemAssembler<Matrix_K1, TrialFunction>& K1_build,
              Matrix_K2& K2_local, SystemAssembler<Matrix_K2, TrialFunctionNeumannSurface>& K2_build)
{

    sp_mat K1;
    K1_build.SetMatrixSize(K); //Settting matrix size here because K1 and K have the same size.
    K1_build.SetMatrixSize(K1);
    K1_build.RunSystemAssembly(K1_local, K1);
    cout<<"Matrix K1 built!\n";
    sp_mat K2;
    K2_build.SetMatrixSize(K2);
    K2_build.RunSystemAssembly(K2_local, K2);
    cout<<"Matrix K2 built!\n";
    K=K1+K2;
    cout<<"Matrix K built!\n";
}

void BuildVector_f(mat& f, SystemAssembler<Matrix_f_Neumann, TrialFunctionNeumannSurface>& f_NeumannBuild,
              Matrix_f_Neumann& fNeumann_local)
{
    f_NeumannBuild.SetVectorSize(f);
    f_NeumannBuild.RunSystemAssemblyVector(fNeumann_local, f);
    cout<<"Vector f built!\n";
}

int main(int argc, char *argv[])
{
    if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./DynamicHeatDiff <.msh Filename> <Dimension>\n";
        return 0;
    }
    double delta_t=1.0,
            thermal_conductivity=50.2,
            conv_ht_coeff=10.0,
            rho=7700.0,
            heat_capacity=510.79,
            AmbientTemperature=25.0,
            BodyTemperature=20.0,
            surf_norm_heat_loss=0.0,
            Applied_Temperature=100.0,
            theta=0.878;
    int vectorLevel=1,
            NeumannBoundaryIndex=0,
            DirichletBoundaryIndex=1,
            BodyIndexForInitialTemp=0;
    std::vector<double> time(2);//A vector of two time variables are take so that 1 takes care of the changes in C, K, etc with time,
                                //the other is for the vector f
    time[0]=0.0;
    time[1]=0.0;
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    libGmshReader::MeshReader Mesh(FileName, Dimension);
    sp_mat C, K;
    std::vector<mat> Sol_of_u(1);
    mat u_ambient;
    mat bar_q_n;
    std::vector<mat> f(2);

    TrialFunction u_dot(Mesh, vectorLevel);
    TestFunctionGalerkin<TrialFunction> v_dot(u_dot);
    FormMultiThread<TrialFunction> a_dot;
    Matrix_C C_local(a_dot, u_dot, v_dot, heat_capacity, rho);
    SystemAssembler<Matrix_C, TrialFunction> C_build(a_dot, u_dot, v_dot);

    BuildMatrix_C(C, C_local, C_build);

    TrialFunction u1(Mesh, vectorLevel);
    TestFunctionGalerkin<TrialFunction> v1(u1);
    FormMultiThread<TrialFunction> a1;
    Matrix_K1 K1_local(a1, u1, v1, thermal_conductivity);
    SystemAssembler<Matrix_K1, TrialFunction> K1_build(a1, u1, v1);
    K1_build.SetVectorSize(Sol_of_u[0]);

    TrialFunctionNeumannSurface u2(u1, NeumannBoundaryIndex);
    TestFunctionGalerkin<TrialFunctionNeumannSurface> v2(u2);
    FormMultiThread<TrialFunctionNeumannSurface> a2;
    Matrix_K2 K2_local(a2, u2, v2, conv_ht_coeff);
    SystemAssembler<Matrix_K2, TrialFunctionNeumannSurface> K2_build(a2, u2, v2);

    BuildMatrix_K(K, K1_local, K1_build, K2_local, K2_build);

    SystemAssembler<Matrix_f_Neumann, TrialFunctionNeumannSurface> f_NeumannBuild(a2, u2, v2);
    std::vector<umat> Local_positions;
    f_NeumannBuild.SetNodePositions_u(Local_positions);
    f_NeumannBuild.SetVectorSize(u_ambient);
    f_NeumannBuild.SetVectorSize(bar_q_n);

    InitialBC_Expression1 TempExpression1(AmbientTemperature, time[1], vectorLevel);
    umat boolAmbientIntialBC={1};
    InitialBC ambient(u2, NeumannBoundaryIndex, boolAmbientIntialBC);
    ambient.SetInitialBCExpression(TempExpression1);
    ambient.ApplyInitialBC(u_ambient);

    InitialBC_Expression1 TempExpression2(BodyTemperature, time[0], vectorLevel);
    InitialBC Body_temp(u1, NeumannBoundaryIndex, boolAmbientIntialBC);
    Body_temp.SetInitialBCExpression(TempExpression2);
    Body_temp.ApplyInitialBC(Sol_of_u[0]);

    InitialBC_Expression2 NormalHeatLossExpression(surf_norm_heat_loss, time[1], vectorLevel);
    InitialBC heatLossNormalSurface(u2, NeumannBoundaryIndex, boolAmbientIntialBC);
    heatLossNormalSurface.SetInitialBCExpression(NormalHeatLossExpression);
    heatLossNormalSurface.ApplyInitialBC(bar_q_n);

    DirichletBC AppliedTemp_Dirichlet(u2, DirichletBoundaryIndex, boolAmbientIntialBC); //reusing boolAmbientInitialBC because Dirichlet BC
    //is applied to the same number of nodes
    DirichletBC_Expression TempExprssn_DiricletBC(Applied_Temperature, vectorLevel);
    AppliedTemp_Dirichlet.SetDirichletBCExpression(TempExprssn_DiricletBC);

    Matrix_f_Neumann f_Neumann_local(a2, u2, v2, bar_q_n, conv_ht_coeff, u_ambient, Local_positions);

    BuildVector_f(f[0], f_NeumannBuild, f_Neumann_local);
    time[1]=(1)*delta_t;
    BuildVector_f(f[1], f_NeumannBuild, f_Neumann_local);
    SingleStepAlgorithm S11_Instance(u1, v1, delta_t);

    GmshWriter GmshWrite(u1, "HtDffssn.msh");
    GmshWrite.viewName="Temperature";

    for (int step=0; step<=100; step++)
    {
        cout<<"Currently at Step "<<step<<"\n";
        sp_mat A;
        mat b;
        GmshWrite.WriteToGmsh(Sol_of_u[0], step, time[0], 0);
        S11_Instance.SetSizeOfMatrix(A);
        S11_Instance.SetSizeOfVector(b);
        S11_Instance.SingleStep_11(A, b, Sol_of_u[0], theta, C, K, f);
        mat alpha;
        AppliedTemp_Dirichlet.ApplyBC_Dynamic(Sol_of_u[0], A, b);
        spsolve(alpha, A, b);
        time[0]=(step+1)*delta_t;
        S11_Instance.update_Sol_of_u(Sol_of_u, alpha);
        time[1]=(step+2)*delta_t;
        mat f_new;
        BuildVector_f(f_new, f_NeumannBuild, f_Neumann_local);
        S11_Instance.update_f(f,f_new);
    }
    return 0;
}
