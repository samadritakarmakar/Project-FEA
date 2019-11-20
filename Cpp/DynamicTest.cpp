#include "ProjectFEA.hpp"

mat get_f(double t)
{
    mat f={pow(t,2)};
    return f;
}

mat get_f2(double t)
{
    mat f={3*exp(t/2.0)};
    return f;
}

int main()
{
    int ch=1;
    cout<<"Enter 1 for solution of dx/dt + 2x =t^2 \n ";
    cout<<"Enter 2 for solution of d^2x/dt^2 + dx/dt = 3 exp(t/2) : ";
    std::cin>>ch;
    if (ch==1)
    {
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
            cout<<"x at t= "<<t+delta_t<<" = "<<x[0](0,0)<<" and alpha = "<<alpha<<"\n";
            double t2=t+delta_t;
            cout<<"Analytical Solution is: x= "<<(pow(t2,2)-t2+0.5)/2.0-13.0/4.0*exp(6.0-2.0*t2)<<"\n";

        }
    }
    else if(ch==2)
    {
        std::vector<mat> x(2);
        x[0].set_size(1,1);
        x[1].set_size(1,1);
        x[0]={4.0};
        x[1]={3};

        std::vector<sp_mat> M_C_K(3);
        for (int i=0; i<M_C_K.size(); i++)
            M_C_K[i].set_size(1,1);
        M_C_K[2](0,0)=1.0;
        M_C_K[1](0,0)=1.0;
        M_C_K[0](0,0)=0.0;

        sp_mat A(1,1); mat b(1,1);
        A.zeros();
        b.zeros();
        double start_t=0.0;
        double max_t=10;
        double delta_t=0.1;
        std::vector<mat> f(3);
        std::vector<double> time(2);
        time[0]=(start_t);
        time[1]=(start_t);
        f[0]=get_f2(time[1]);
        time[1]+=delta_t;
        f[1]=get_f2(time[1]);
        time[1]+=delta_t;
        f[2]=get_f2(time[1]);
        SingleStepAlgorithm TestAlgo(delta_t);
        std::vector<double> theta(2);
        theta[0]=0.5;
        theta[1]=0.5;
        for (int step=0; step<=max_t/delta_t; step++)
        {
            TestAlgo.SingleStep_22(A, b, x, theta, M_C_K, f);
            mat alpha;
            spsolve(alpha,A,b);
            TestAlgo.update_Sol_of_u(x, alpha);
            time[1]+=delta_t;
            mat f_new=get_f2(time[1]);
            TestAlgo.update_f(f, f_new);
            cout<<"x at t= "<<time[0]+delta_t<<" = "<<x[0](0,0)<<" and alpha = "<<alpha;
            double t2=time[0]+delta_t;
            cout<<"Analytical Sol at t = "<<t2<<" = "<<-exp(-t2)+1.0+4.0*exp(t2/2.0)<<"\n\n";
            time[0]+=delta_t;
        }
    }
    return 0;
}
