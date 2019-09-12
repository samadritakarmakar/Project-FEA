#include "ProjectFEA.hpp"
#include <math.h>

void GetStabilizationParameters(const double& sizeElement, const vec& vel, const double& nu, const double& sigma, double& tau)
{
    double normVel=norm(vel);
    double sizeElement2=sizeElement;
    double Pe=normVel*sizeElement2/(2.0*nu);
    double pe1=vel(0)*sizeElement2/(2.0*nu);
    double pe2=vel(1)*sizeElement2/(2.0*nu);
    Pe=(pe1+pe2)/2;
    double beta=(1.0/(1.0/std::tanh(Pe)-1.0/Pe));
    //tau=(sizeElement/(2.0*normVel))*1.0/(1.0+1.0/Pe+sizeElement*sigma/(2.0*normVel));
    //tau.set_size(vel.n_rows, vel.n_rows);
    tau=beta*sizeElement2/(2.0*normVel);
    tau=1.0;
}

class FindSize : public LocalIntegrator<TrialFunction>
{
public:
    FindSize(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v):
        LocalIntegrator<TrialFunction> (a, u, v){}
    double scalar_integration(FormMultiThread<TrialFunction>& a, TrialFunction& u, int thread)
    {
        return a[thread].dX(u);
    }
};

class DiffusionLHS : public LocalIntegrator<TrialFunction>
{
public:
    DiffusionLHS(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v,
                 vec& h, vec& vel, double& nu, double& sigma):
        LocalIntegrator<TrialFunction> (a, u, v),
        h(h), nu(nu), sigma(sigma), vel(vel){}

    sp_mat weak_form(FormMultiThread<TrialFunction>& a, TrialFunction& u,
                         TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        sp_mat LHS;
        sp_mat R_u_LHS;
        sp_mat P_v_LHS;
        double hTemp; //Pe;
        double tau;
        /// The size of the Element. Volume for 3D, Area for 2D.
        hTemp=(h(a[thread].ElementNumber));
        GetStabilizationParameters(hTemp, vel, nu, sigma, tau);
        ///VariableGeneric<double>R_u_RHS(a.GetNumOfThreads());
        /// Represents (a.grad(u)+sigma*u) -------- See equation 2.55 in Finite Element Methods
        /// for Flow Problems. Jean Donea and Antonio Huerta.
        R_u_LHS=a[thread].vctr_dot_grad_u(vel,u)+sigma*a[thread].u(u);
        //GLS Stabilization for Linear Elements
        P_v_LHS=a[thread].vctr_dot_grad_v(vel,v)+sigma*a[thread].v(v);
        //Represents nu*dot(grad(v),grad(u))+v*dot(vel,grad(u))+sigma*u)
        LHS=nu*a[thread].dot(a[thread].grad(v),a[thread].grad(u))+
                a[thread].v(v)*(a[thread].vctr_dot_grad_u(vel,u)+sigma*a[thread].u(u));

        return (LHS+P_v_LHS*tau*R_u_LHS)*a[thread].dX(u);
    }
    vec &h, &vel;
    double &nu;
    double &sigma;
};
class DiffusionRHS : public LocalIntegrator<TrialFunction>
{
public:
    DiffusionRHS(FormMultiThread<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v,
                 vec& h, vec& vel, double& nu, double& sigma, double& source):
        h(h), vel(vel), nu(nu), sigma(sigma), source(source),
        LocalIntegrator<TrialFunction> (a, u, v){}

    mat weak_form_vector(FormMultiThread<TrialFunction>& a, TrialFunction& u,
                             TestFunctionGalerkin<TrialFunction>& v, int thread)
    {
        //double Pe;
        double tau;
        double hTemp=(h(a[thread].ElementNumber));
        GetStabilizationParameters(hTemp, vel, nu, sigma, tau);
        double R_u_RHS=source;
        sp_mat P_v_RHS=a[thread].vctr_dot_grad_v(vel,v)+sigma*a[thread].v(v);;
        sp_mat RHS=a[thread].v(v)*source;
        return mat(P_v_RHS*tau*R_u_RHS+RHS)*a[thread].dX(u);
    }
    vec &h, &vel;
    double &nu;
    double &sigma, &source;
};

class NormalFlux : public LocalIntegrator<TrialFunctionNeumannSurface>
{
public:
    NormalFlux(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
               TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, double& NormFluxVel):
        LocalIntegrator<TrialFunctionNeumannSurface> (a, u, v), h_T(NormFluxVel) {}

    mat weak_form_vector(FormMultiThread<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u,
                         TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, int thread)
    {
        return h_T*mat(a[thread].v(v))*a[thread].dS(u);
    }
    double h_T;
};


class NormalFluxLine : public LocalIntegrator<TrialFunctionNeumannLine>
{
public:
    NormalFluxLine(FormMultiThread<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
               TestFunctionGalerkin<TrialFunctionNeumannLine>& v, double& NormFluxVel):
        LocalIntegrator<TrialFunctionNeumannLine> (a, u, v), h_T(NormFluxVel) {}

    mat weak_form_vector(FormMultiThread<TrialFunctionNeumannLine>& a, TrialFunctionNeumannLine& u,
                         TestFunctionGalerkin<TrialFunctionNeumannLine>& v, int thread)
    {
        return h_T*mat(a[thread].v(v))*a[thread].dL(u);
    }
    double h_T;
};

///This class over here through its overloaded virtual function declares the values of Dirichlet Nodes.
/// The virtual function 'Eval' is evaluated at each node to find the value of Dirichlet Condtion at that node.
class DirichletExprssn : public Expression
{
public:
    DirichletExprssn (int vectorLevel): Expression (vectorLevel)
    {
    }

    vec Eval(vec& x)
    {
        vec value={0};
        return value;
    }
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

    //Parameters of Advection-Diffusion-Reaction equation.
    //Diffusion coefficent of Carbon Monoxoide in air
    double nu=0.208e-4;
    double sigma=0.0, source=0.0, NrmlFlux=-5.0;
    vec vel;
    if(Dimension==3)
    {
        vel<<0.0<<endr<<4.0<<endr<<1.0;
    }
    else if (Dimension==2)
    {
        vel<<1.0<<endr<<4.0;
    }
    TestFunctionGalerkin<TrialFunction> v(u);
    FormMultiThread<TrialFunction> a;
    vec h;
    SystemAssembler<FindSize,TrialFunction> SclrIntrgtn(a, u ,v);
    FindSize SizeOfElement(a,u,v);
    SclrIntrgtn.RunScalarIntegration(SizeOfElement, h);
    h=abs(h);
    cout<<"Total Volume = "<<sum(h)<<"\n";

    SystemAssembler<DiffusionLHS, TrialFunction> Dffsn_LHS_MtrxSystm(a, u, v);
    DiffusionLHS Dffsn_LHS(a, u, v, h, vel, nu, sigma);
    sp_mat A;
    Dffsn_LHS_MtrxSystm.SetMatrixSize(A);
    Dffsn_LHS_MtrxSystm.RunSystemAssembly(Dffsn_LHS, A);

    SystemAssembler<DiffusionRHS, TrialFunction> Dffsn_RHS_MtrxSystm(a, u, v);
    DiffusionRHS Dffsn_RHS(a, u, v, h, vel, nu, sigma, source);
    mat b;
    Dffsn_RHS_MtrxSystm.SetVectorSize(b);
    Dffsn_RHS_MtrxSystm.RunSystemAssemblyVector(Dffsn_RHS, b);
   // cout<<"b= after Diffusion RHS"<<b;

    /*TrialFunctionNeumannSurface u1(u,0);
    TestFunctionGalerkin<TrialFunctionNeumannSurface> v1(u1);
    FormMultiThread<TrialFunctionNeumannSurface> a1;
    NormalFlux NrmalFlx(a1,u1,v1,NrmlFlux);
    SystemAssembler<NormalFlux, TrialFunctionNeumannSurface> NrmlFluxVctrSystm(a1,u1,v1);
    NrmlFluxVctrSystm.RunSystemAssemblyVector(NrmalFlx, b);*/
    //cout<<"b= after NormalFlux"<<b;

    TrialFunctionNeumannLine u2(u,0);
        TestFunctionGalerkin<TrialFunctionNeumannLine> v2(u2);
        FormMultiThread<TrialFunctionNeumannLine> a2;
        NormalFluxLine NrmalFlx2(a2,u2,v2,NrmlFlux);
        SystemAssembler<NormalFluxLine, TrialFunctionNeumannLine> NrmlFluxVctrSystm(a2,u2,v2);
        NrmlFluxVctrSystm.RunSystemAssemblyVector(NrmalFlx2, b);

    umat boolDirichletNodes={1};
    DirichletBC DrcltBC(u2,1,boolDirichletNodes);
    DirichletExprssn Exprssn(1);
    DrcltBC.SetDirichletBCExpression(Exprssn);
    DrcltBC.ApplyBC(A,b);

    //cout<<"A="<<mat(A);
    //cout<<"b= "<<b;

    vec x;


    x=spsolve(A, b);
    //cout<<"x="<<x;
    GmshWriter WriteDensityField(u, "DnstyFld.pos");
    WriteDensityField.WriteToGmsh(x);
}
