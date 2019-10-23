#include "ProjectFEA.hpp"

mat get_f(double t)
{
    mat f={pow(t,2)};
    return f;
}

int main(int argc, char *argv[])
{
    if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./DiffusionEx <.msh Filename> <Dimension>\n";
        return 0;
    }
    std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    libGmshReader::MeshReader Mesh(FileName, Dimension);
    TrialFunction dummy(Mesh, 1);
    TestFunctionGalerkin<TrialFunction> dummy_v(dummy);
    mat x(1,1);
    sp_mat A(1,1);
    mat b(1,1);
    sp_mat C(1,1);
    C(0,0)=1.0;
    sp_mat K(1,1);
    K.zeros();
    K(0,0)=2.0;
    mat x0={0.0};
    x=x0;
    double start_t=3.0;
    double max_t=1000;
    double delta_t=0.1;
    SingleStepAlgorithm TestAlgo(dummy, dummy_v, delta_t);
    TestAlgo.SetSingleStep_theta_To(0.5);
    mat fn0=get_f(0.0);
    for (double t=start_t; t<max_t; t+=delta_t)
    {
        mat fn1=get_f(t+delta_t);
        TestAlgo.SingleStep_11(A, b, x, 0.5, C, K, fn0, fn1);
        mat alpha;
        spsolve(alpha, A, b);
        x=x+delta_t*alpha;
        fn0=fn1;
        cout<<"x at t= "<<t+delta_t<<" = "<<x<<" and alpha = "<<alpha<<"\n";
    }
}
