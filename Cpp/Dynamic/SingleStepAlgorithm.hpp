#ifndef SINGLESTEPALGORITHM_HPP
#define SINGLESTEPALGORITHM_HPP
#include "TrialFunction.hpp"
#include "TestFunctionGalerkin.hpp"
#include "ProjectFEA_Math.hpp"
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

    void SingleStep_11(sp_mat& A, mat& b, mat& Sol_of_u_n0, double theta, sp_mat& C_Of_dot_u, sp_mat& K_Of_u,  mat& f_n0, mat& f_n1)
    {
        std::vector<sp_mat> Mat_variables(2);
        Mat_variables[1].set_size(C_Of_dot_u.n_rows,C_Of_dot_u.n_cols);
        Mat_variables[1]=C_Of_dot_u;
        Mat_variables[0].set_size(K_Of_u.n_rows,K_Of_u.n_cols);
        Mat_variables[0]=K_Of_u;
        std::vector<mat> f(2);
        f[0].set_size(f_n0.n_rows, f_n0.n_cols);
        f[0]=f_n0;
        f[1].set_size(f_n1.n_rows, f_n1.n_cols);
        f[1]=f_n1;
        std::vector<mat> Sols_of_u(1);
        Sols_of_u[0]=Sol_of_u_n0;
        //Sols_of_u[1]=Sol_of_u_n1;
        General_SSpj(A, b, Mat_variables, Sols_of_u, f);
    }

    void SetSingleStep_theta_To(double theta_ss)
    {
        theta.push_back(theta_ss);
        theta_f.push_back(theta_ss);
    }

    void Set_theta_To(std::vector<double>& theta_vector)
    {
        theta=theta_vector;
    }

    void Set_theta_f_To(std::vector<double>& theta_f_vector)
    {
        theta_f=theta_f_vector;
    }


    /// Implementation of Single Step time algorithm as per O.C. Zienkiewicz book.
    /// A, b are the values returned, Mat_variables is a vector with index 0 as K, 1 as C and 2 as M; or 0 as K, 1 as C
    /// size_u_mean refers to the degree of u_mean. Should not be greater than 3.
    void General_SSpj(sp_mat& A, mat& b, std::vector<sp_mat>& Mat_variables, std::vector<mat>& Sols_of_u, std::vector<mat>& f)
    {
        std::vector<mat>& u=Sols_of_u;
        int size_u_mean = u.size();
        int p=theta.size();
        if(p!=u.size())
        {
            cout<<"The vector size of theta has to be equal to the vector size of u\n";
            throw;
        }
        /*if(size_u_mean>3 || size_u_mean<1)
        {
            cout<<"size_u_mean cannot be less than 1 and greater than 3\n";
            throw;
        }*/
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
        for(int dot=1; dot<theta.size(); dot++)
        {
            A=A+theta[dot]*pow(delta_t, dot)*Mat_variables[k]/double(factorial(dot));
            k--;
        }
        mat f_mean;
        GetMean_f(f, f_mean);
        b=-b+f_mean;
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
