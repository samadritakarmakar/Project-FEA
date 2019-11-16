#ifndef SINGLESTEPALGORITHM_HPP
#define SINGLESTEPALGORITHM_HPP
#include "TrialFunction.hpp"
#include "TestFunctionGalerkin.hpp"
#include "ProjectFEA_Math.hpp"
#include "DirichletBC.hpp"
//Author of code Samadrita Karmakar.
/// Refer The Finite Element Method: Its Basis and Fundamentals; Sixth edition
/// by O.C. Zienkiewicz, R.L. Taylor, J.Z. Zhu for this implementation



class SingleStepAlgorithm
{
public:
    SingleStepAlgorithm(TrialFunction& u, TestFunctionGalerkin<TrialFunction>& v, double delta_t): delta_t(delta_t)
    {
        for (int ElmntType=0; ElmntType<u.NoOfElementTypes; ElmntType++)
        {
            no_of_rows+=u.numOfNodes*u.vectorLvl;
            no_of_cols+=v.numOfNodes*v.vectorLvl;
        }
    }

    SingleStepAlgorithm(double delta_t): delta_t(delta_t)
    {
    }

    /// This function uses a specific case of General_SSpj, where p=1 & j=1.
    /// A and b are the returned matrices. Sol_of_u_n0 is the current value of the solution.
    /// theta >0.5 gives a stable algorithm.
    /// f is a vector of type mat. f[0] is f_n, f[1] is f_n+1, f[2] is f_n+2, etc.
    void SingleStep_11(sp_mat& A, mat& b, mat& Sol_of_u_n0, double theta, sp_mat& C_Of_dot_u, sp_mat& K_Of_u,  std::vector<mat>& f,
                       DirichletBC & applied_Dirichlet)
    {
        std::vector<sp_mat> Mat_variables(2);
        Mat_variables[1].set_size(C_Of_dot_u.n_rows,C_Of_dot_u.n_cols);
        Mat_variables[1]=C_Of_dot_u;
        Mat_variables[0].set_size(K_Of_u.n_rows,K_Of_u.n_cols);
        Mat_variables[0]=K_Of_u;
        SetSingleStep_11_theta_To(theta);
        std::vector<mat> Sols_of_u(1);
        Sols_of_u[0]=Sol_of_u_n0;
        General_SSpj(A, b, Mat_variables, Sols_of_u, f, applied_Dirichlet);
    }

    void SingleStep_11(sp_mat& A, mat& b, mat& Sol_of_u_n0, double theta, sp_mat& C_Of_dot_u, sp_mat& K_Of_u,  std::vector<mat>& f)
    {
        std::vector<sp_mat> Mat_variables(2);
        Mat_variables[1].set_size(C_Of_dot_u.n_rows,C_Of_dot_u.n_cols);
        Mat_variables[1]=C_Of_dot_u;
        Mat_variables[0].set_size(K_Of_u.n_rows,K_Of_u.n_cols);
        Mat_variables[0]=K_Of_u;
        SetSingleStep_11_theta_To(theta);
        std::vector<mat> Sols_of_u(1);
        Sols_of_u[0]=Sol_of_u_n0;
        General_SSpj(A, b, Mat_variables, Sols_of_u, f);
    }

    void SetSingleStep_11_theta_To(double theta_ss)
    {
        theta=std::vector<double>(1);
        theta_f=std::vector<double>(1);
        theta[0]=theta_ss;
        theta_f[0]=theta_ss;
    }

    void Set_theta_To(std::vector<double>& theta_vector)
    {
        theta=theta_vector;
    }

    void Set_theta_f_To(std::vector<double>& theta_f_vector)
    {
        theta_f=theta_f_vector;
    }

    void update_f(std::vector<mat>& current_f, mat& f_new)
    {
        int i=0;
        for(i=0; i<current_f.size()-1; i++)
        {
            current_f[i]=current_f[i+1];
        }
        current_f[i]=f_new;
    }

    void update_Sol_of_u(std::vector<mat>& Sols_of_u, mat& alpha)
    {
        for (int dot=0; dot<Sols_of_u.size(); dot++)
        {
            int final_n=1;
            for(int taylor_n=dot+1; taylor_n<Sols_of_u.size(); taylor_n++)
            {
                Sols_of_u[dot]=Sols_of_u[dot]+Sols_of_u[taylor_n]*pow(delta_t,taylor_n)/factorial(taylor_n);
                final_n++;
            }
            Sols_of_u[dot]=Sols_of_u[dot]+alpha*pow(delta_t, final_n)/factorial(final_n);
        }
    }

    /// Implementation of Single Step time algorithm as per O.C. Zienkiewicz book.
    /// A, b are the values returned.
    /// Mat_variables is a vector with index 0 as K, 1 as C and 2 as M; or 0 as K, 1 as C.
    /// Sols_of_u is a vector of the solution matrix of u with u[0]=u, u[1]=\dot{u}, u[2]=\ddot{u}, etc.
    /// f is a vector of type mat. f[0] is f_n, f[1] is f_n+1, f[2] is f_n+2, etc.
    void General_SSpj(sp_mat& A, mat& b, std::vector<sp_mat>& Mat_variables, std::vector<mat>& Sols_of_u, std::vector<mat>& f,
                      DirichletBC & applied_Dirichlet)
    {
        /*std::vector<mat>& u=Sols_of_u;
        int size_u_mean = u.size();
        int p=theta.size();
        if(p!=u.size())
        {
            cout<<"The vector size of theta has to be equal to the vector size of u\n";
            throw;
        }
        std::vector<mat> u_mean(size_u_mean);
        //SetSizeOfMatrix(A);
        //SetSizeOfVector(b);
        A.zeros();
        b.zeros();
        for (int dot=0; dot<=size_u_mean-1; dot++)
        {
            u_mean[dot].set_size(u[0].n_rows,u[0].n_cols);
            u_mean[dot]=u[dot];
            for (int q=dot; q<p-1; q++)
            {
                u_mean[dot]=u_mean[dot]+theta[q]*pow(delta_t,q+1)*u[q+1]/double(factorial(q+1-dot));
            }
            b=b+Mat_variables[dot]*u_mean[dot];
        }
        A=A+Mat_variables[theta.size()];
        int k=theta.size();
        for(int dot=1; dot<=theta.size(); dot++)
        {
            k--;
            //cout<<"Mat_variables[k] at "<<k<<" = "<<Mat_variables[k]<<"\n" ;
            A=A+theta[dot-1]*pow(delta_t, dot)*Mat_variables[k]/double(factorial(dot));
        }
        mat f_mean;
        GetMean_f(f, f_mean);
        //cout<<"f_mean = "<<f_mean<<"\n";
        b=-b+f_mean;*/
        Build_b(b, Mat_variables, Sols_of_u, f);
        Build_A(A, Mat_variables, Sols_of_u);
        ApplyDirichletBC_to_M_C_K(Mat_variables, applied_Dirichlet);
    }

    void General_SSpj(sp_mat& A, mat& b, std::vector<sp_mat>& Mat_variables, std::vector<mat>& Sols_of_u, std::vector<mat>& f)
    {
        Build_b(b, Mat_variables, Sols_of_u, f);
        Build_A(A, Mat_variables, Sols_of_u);
    }

    void GetMean_f(std::vector<mat>& f, mat& f_mean)
    {
        if(f.size()!=theta_f.size()+1)
        {
            cout<<"Size of std::vector f must be 1 greater than the size of theta_f\n";
            throw;
        }
        //SetSizeOfVector(f_mean);
        f_mean.set_size(f[0].n_rows, f[0].n_cols);
        std::vector<mat> df_by_dt;
        Calculate_df_by_dt(f, df_by_dt);
        f_mean=f[0];
        for(int j=0; j<theta_f.size(); j++)
        {
            f_mean=f_mean+theta_f[j]*pow(delta_t,j+1)*df_by_dt[j];
            //cout<<"df_by_dt = "<<df_by_dt[j]<<" at j= "<<j<<"\n";
        }
    }

    void SetSizeOfMatrix(sp_mat& A)
    {
        A.set_size(no_of_rows, no_of_cols);
    }
    void SetSizeOfVector(mat& b)
    {
        b.zeros(no_of_rows, 1);
    }
private:
    int no_of_rows=0;
    int no_of_cols=0;
    double delta_t;
    std::vector<double> theta;
    std::vector<double> theta_f;

    void ApplyDirichletBC_to_M_C_K(std::vector<sp_mat>& Mat_variables, DirichletBC & applied_Dirichlet)
    {
        for(int i=0; i<Mat_variables.size(); i++)
        {
            applied_Dirichlet.Set_Dynamic_BC_on_M_K_C(Mat_variables[i]);
        }
    }

    void Build_b(mat &b, std::vector<sp_mat>& Mat_variables, std::vector<mat>& Sols_of_u, std::vector<mat>& f)
    {
        std::vector<mat>& u=Sols_of_u;
        int size_u_mean = u.size();
        int p=theta.size();
        if(p!=u.size())
        {
            cout<<"The vector size of theta has to be equal to the vector size of u\n";
            throw;
        }
        std::vector<mat> u_mean(size_u_mean);
        //SetSizeOfMatrix(A);
        //SetSizeOfVector(b);
        for (int dot=0; dot<=size_u_mean-1; dot++)
        {
            u_mean[dot].set_size(u[0].n_rows,u[0].n_cols);
            u_mean[dot]=u[dot];
            for (int q=dot; q<p-1; q++)
            {
                u_mean[dot]=u_mean[dot]+theta[q]*pow(delta_t,q+1)*u[q+1]/double(factorial(q+1-dot));
            }
            b=b+Mat_variables[dot]*u_mean[dot];
        }
        mat f_mean;
        GetMean_f(f, f_mean);
        //cout<<"f_mean = "<<f_mean<<"\n";
        b=-b+f_mean;
    }

    void Build_A(sp_mat &A, std::vector<sp_mat>& Mat_variables, std::vector<mat>& Sols_of_u)
    {
        A=A+Mat_variables[theta.size()];
        int k=theta.size();
        for(int dot=1; dot<=theta.size(); dot++)
        {
            k--;
            //cout<<"Mat_variables[k] at "<<k<<" = "<<Mat_variables[k]<<"\n" ;
            A=A+theta[dot-1]*pow(delta_t, dot)*Mat_variables[k]/double(factorial(dot));
        }
    }

    void Calculate_df_by_dt(std::vector<mat>& f, std::vector<mat>& df_by_dt)
    {
        mat Delta_t_Matrix;
        df_by_dt=std::vector<mat> (theta_f.size());
        BuildDelta_t_Matrix(Delta_t_Matrix);
        mat diff_f(theta_f.size(),f[0].n_rows);
        for (int j=0; j<theta_f.size(); j++)
        {
            diff_f.row(j)=f[j+1].t()-f[0].t();
        }
        mat df_dt;
        solve(df_dt, Delta_t_Matrix, diff_f);
        for(int j=0; j<df_dt.n_rows; j++)
        {
            df_by_dt[j].set_size(df_dt.n_cols,1);
            df_by_dt[j] =df_dt.row(j).t();
        }
    }

    void BuildDelta_t_Matrix(mat& Delta_t_Matrix)
    {
        Delta_t_Matrix.set_size(pow(theta_f.size(),2));
        for (int i=0; i<theta_f.size(); i++)
        {
            Delta_t_Matrix(i,0)=(i+1)*delta_t;
        }
        for(int j=1; j<theta_f.size(); j++)
        {
            Delta_t_Matrix.col(j)=pow(Delta_t_Matrix.col(0),j+1)/double(factorial(j+1));
        }
    }
};

#endif // SINGLESTEPALGORITHM_HPP
