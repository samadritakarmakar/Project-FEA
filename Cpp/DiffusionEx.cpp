#include "ProjectFEA.hpp"
class FindSize : public LocalIntegrator<TrialFunction>
{
public:
    FindSize(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator<TrialFunction> (a, u, v){}
    double scalar_integration(FormMultiThread<TrialFunction>& a, TrialFunction& u, int thread)
    {
       // if (a.ElementNumber==134)
     //  {
            //cout<<"a[thread].ElementNumber= "<<a[thread].ElementNumber<<" ";
     //       cout<<"a[thread].dX(u)" <<a.dX(u)<<" ";
      // }
        return a[thread].dX(u);
    }
};

class DiffusionLHS : public LocalIntegrator<TrialFunction>
{
public:
    DiffusionLHS(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator<TrialFunction> (a, u, v){}
    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u,
                         TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        //Diffusion coefficent of Carbon Monoxoide in air
        double nu=0.208e-4;
        double sigma=10.0;
        vel=10;
        //double source=20.0;
        //GLS Stabilization
        VariableGeneric<sp_mat> R_u_LHS(a.GetNumOfThreads());
        VariableGeneric<sp_mat> P_w_LHS(a.GetNumOfThreads());
        VariableGeneric<double> tau(a.GetNumOfThreads());
        //The size of the Element. Volume for

        //VariableGeneric<double>R_u_RHS(a.GetNumOfThreads());
        // Represents (a.grad(u)-sigma*u) -------- See equation 2.61 in Finite Element Methods
        //for Flow Problems. Jean Donea and Antonio Huerta.
        R_u_LHS[thread]=a[thread].dot(vel, a[thread].grad(u))+sigma*a[thread].u(u);
        P_w_LHS[thread]=a[thread].grad(v);
        //R_u_RHS[thread]=source;
        //Represents v*(vel.grad(u))dX
        return a[thread].v(v)*a[thread].dot(vel,a[thread].grad(u))*a[thread].dX(u);
    }
    vec vel;
};
class DiffusionRHS : public LocalIntegrator<TrialFunction>
{
public:
    DiffusionRHS(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator<TrialFunction> (a, u, v){}
    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u,
                         TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        //Diffusion coefficent of Carbon Monoxoide in air
        double nu=0.208e-4;
        double sigma=10.0;
        vel=10;
        //double source=20.0;
        //GLS Stabilization
        VariableGeneric<sp_mat> R_u_LHS(a.GetNumOfThreads());
        //R_u_RHS[thread]=source;
        return a[thread].v(v)*a[thread].dot(vel,a[thread].grad(u))*a[thread].dX(u);
    }
    vec vel;
};

int main(int argc, char *argv[])
{
    if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./DiffusionEx <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    int vectorLevel=1;
    libGmshReader::MeshReader Mesh(FileName, Dimension);
    libGmshReader::MeshReader Mesh_order_1(Mesh, Dimension, 1);
    TrialFunction u(Mesh_order_1, vectorLevel);
    //TrialFunction u(Mesh, vectorLevel);
    TestFunctionGalerkin<TrialFunction> v(u);
    FormMultiThread<TrialFunction> a;
    vec h;
    SystemAssembler<FindSize,TrialFunction> SclrIntrgtn(a, u ,v);
    FindSize SizeOfElement(a,u,v);
    SclrIntrgtn.RunScalarIntegration(SizeOfElement, h);
    cout<<"Total Volume = "<<sum(h)<<"\n";
}
