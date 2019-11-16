#include "ProjectFEA.hpp"

mat get_f(double t)
{
    mat f={pow(t,2)};
    return f;
}

int main(int argc, char *argv[])
{
    /*if(argc==1 || argc>3)
    {
        std::cout<<"Usage: ./TestDynamic <.msh Filename> <Dimension>\n";
        return 0;
    }*/
    /*std::string FileName(argv[1]);
    int Dimension=*argv[2]-'0';
    libGmshReader::MeshReader Mesh(FileName, Dimension);
    TrialFunction dummy(Mesh, 1);
    TestFunctionGalerkin<TrialFunction> dummy_v(dummy);*/
    std::vector<mat> x(1);
    x[0].set_size(1,1);
    sp_mat A(1,1);
    mat b(1,1);
    sp_mat C(1,1);
    C(0,0)=1.0;
    sp_mat K(1,1);
    K.zeros();
    K(0,0)=2.0;
    mat x0={0.0};
    x[0]=x0;
    double start_t=3.0;
    double max_t=10;
    double delta_t=0.1;
    SingleStepAlgorithm TestAlgo(delta_t);
    //SingleStepAlgorithm TestAlgo(dummy, dummy_v, delta_t);
    TestAlgo.SetSingleStep_11_theta_To(0.5);
    std::vector<mat> f(2);
    f[0]=get_f(start_t);
    f[1]=get_f(start_t+delta_t);
    for (double t=start_t; t<max_t; t+=delta_t)
    {
        //mat fn1=get_f(t+delta_t);
        TestAlgo.SingleStep_11(A, b, x[0], 0.5, C, K, f);
        mat alpha;
        cout<<"A = "<<mat(A)<<"\n";
        spsolve(alpha, A, b);
        TestAlgo.update_Sol_of_u(x, alpha);
        mat f_new=get_f(t+delta_t+delta_t); //Gets ready to put new value of f_n+2.
        TestAlgo.update_f(f, f_new);
        //x=x+delta_t*alpha;
        //fn0=fn1;
        cout<<"x at t= "<<t+delta_t<<" = "<<x[0]<<" and alpha = "<<alpha<<"\n";
    }
}
