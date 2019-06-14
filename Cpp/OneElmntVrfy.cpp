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
        cout<<"dN_by_dEps at "<<a.GaussPntr<<" = \n"<<u.dN_by_dEps[a.ElementType][a.GaussPntr]<<"\n";
        cout<<"grad(u) at gauss pt "<<a.GaussPntr<<"=\n "<<mat(a.grad(u));
        cout<<"grad(v):grad(u) at gauss pt "<<a.GaussPntr<<"=\n "<<mat(a.inner(a.grad(v),a.grad(u)))<<"\n";
        cout<<"dX at gauss pt "<<a.GaussPntr<<"= "<<a.dX(u)<<"\n";
        mat F;
        u.Get_F(a.ElementType,a.ElementNumber,a.GaussPntr,F);
        cout<<"F=\n"<<F<<"\n";
        return a.inner(a.grad(v),a.grad(u))*a.dX(u);
    }
};

int main()
{
    libGmshReader::MeshReader Mesh("OneElmntGeom/Triangle.msh",2);
    TrialFunction u(Mesh, 1);
    TestFunctionGalerkin<TrialFunction> v(u);
    Form<TrialFunction> a(u);
    sp_mat A;
    new_LocalIntegrator intgrt(a,u,v);
    SystemAssembler<new_LocalIntegrator, TrialFunction> systmAssmly(a,u,v);
    systmAssmly.SetMatrixSize(A);
    systmAssmly.RunSystemAssembly(intgrt,A);
    cout<<"A= \n"<<mat(A);

    return 0;
}
