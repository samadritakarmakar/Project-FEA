#include "ProjectFEA.hpp"
class Matrix_K: public LocalIntegrator<TrialFunction>
{
public:
    Matrix_K(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, double thermal_conductivity):
        LocalIntegrator (a,u,v), k(thermal_conductivity)
    {
    }
    sp_mat weak_form(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v)
    {
        return k*a.inner(a.grad(v),a.grad(u))*a.dX(u);
    }
    double k;
};

class Matrix_C: public LocalIntegrator<TrialFunction>
{
public:
    Matrix_C(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, double heatCapacity, double density):
        LocalIntegrator (a,u,v), c(heatCapacity), rho(density)
    {
    }
    sp_mat weak_form(Form<TrialFunction>& a, TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v)
    {
        return rho*c*a.v(v)*a.u(u)*a.dX(u);
    }
    double c, rho;
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
        return mat(a.v(v)*s*a.dS(u));
    }
    double& s;
};

class Matrix_f_Neumann: public LocalIntegrator<TrialFunctionNeumannSurface>
{
public:
    Matrix_f_Neumann(Form<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u, TestFunctionGalerkin<TrialFunctionNeumannSurface>& v, double &heat_normal):
        LocalIntegrator (a,u,v), q_bar(heat_normal)
    {
    }
    mat weak_form_vector(Form<TrialFunctionNeumannSurface>& a, TrialFunctionNeumannSurface& u, TestFunctionGalerkin<TrialFunctionNeumannSurface>& v)
    {
        return mat(a.v(v)*q_bar*a.dS(u));
    }
    double& q_bar;
};

