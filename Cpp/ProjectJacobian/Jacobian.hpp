#ifndef JACOBIAN_HPP
#define JACOBIAN_HPP
#include <armadillo>
using namespace arma;

/// Stores x_h in variable 'vec x_h' and returns a matrix of the form
/// [x+x_h, y, z; x, y+y_h, z; x, y, z+z_h]
inline mat Get_xPlus_h(vec x, vec& x_h);

template <class classIntnc, class OtherData>
mat Jacobian(classIntnc *Instnc, mat (classIntnc::*f)(mat , OtherData), mat x1, OtherData data);

template<class OtherData>
mat Jacobian(mat (*f)(mat , OtherData), mat x1, OtherData data) //pointer of the function; and the column of vectors
{
    int i;
    vec x=vectorise(x1); //Ensures a column vector
    vec x_h;
    mat xPlus_h=Get_xPlus_h(x, x_h); //Matrix of xPlus_h and values of x_h
    const vec fconst= vectorise((*f)(x, data)); //values for f(x)
    vec fx1=vectorise((*f)(xPlus_h.row(0), data))-fconst; //for the 1st row; (f(x+x_h) - f(x))
    mat df_dx=zeros(fx1.n_rows, x.n_rows);
    for (i=0; i<x.n_rows-1; i++)
    {
        df_dx.col(i)=fx1/x_h(i); //(f(x+x_h)(i) - f(x))/x_h(i)
        fx1=vectorise((*f)(xPlus_h.row(i+1), data))-fconst;
    }
    df_dx.col(i)=fx1/x_h(i); //for the last column

    return df_dx;
}

template <class classIntnc, class OtherData>
mat Jacobian(classIntnc *Instnc, mat (classIntnc::*f)(mat , OtherData), mat x1, OtherData data) //pointer of the function; and the column of vectors
{
    int i;
    vec x=vectorise(x1); //Ensures a column vector
    vec x_h;
    mat xPlus_h=Get_xPlus_h(x, x_h); //Matrix of xPlus_h and values of x_h
    const vec fconst= vectorise((Instnc->*f)(x, data)); //values for f(x)
    vec fx1=vectorise((Instnc->*f)(xPlus_h.row(0), data))-fconst; //for the 1st row; (f(x+x_h) - f(x))
    mat df_dx=zeros(fx1.n_rows, x.n_rows);
    for (i=0; i<x.n_rows-1; i++)
    {
        df_dx.col(i)=fx1/x_h(i); //(f(x+x_h)(i) - f(x))/x_h(i)
        fx1=vectorise((Instnc->*f)(xPlus_h.row(i+1), data))-fconst;
    }
    df_dx.col(i)=fx1/x_h(i); //for the last column

    return df_dx;
}

template <class Mat1, class OtherData>
Mat1 JacobianMatrix(Mat1 (*f)(mat , OtherData), mat x1, OtherData data)
{
    vec x=vectorise(x1); //Ensures a column vector
    vec x_h;
    mat xPlus_h=Get_xPlus_h(x, x_h); //Matrix of xPlus_h and values of x_h
    const Mat1 fconst=(*f)(x, data); //Change to Mat1
    int NoOfRows=fconst.n_rows;
    int NoOfColumns=fconst.n_cols;
    Mat1 df_dx(NoOfRows, NoOfColumns*x.n_rows); //Change to Mat1
    for (int i=0; i<x.n_rows; i++)
    {
        Mat1 fx1=(*f)(xPlus_h.row(i), data)-fconst;
        Mat1 dfx1_dx=fx1/x_h(i);
        for (int j=0; j<NoOfColumns; j++)
        {
            df_dx.col(NoOfColumns*j+i)=dfx1_dx.col(j);
        }
    }
    return df_dx;
}

template <class classIntnc, class Mat1, class OtherData>
Mat1 JacobianMatrix(classIntnc *Instnc, Mat1 (classIntnc::*f)(mat , OtherData), mat x1, OtherData data)
{
    vec x=vectorise(x1); //Ensures a column vector
    vec x_h;
    mat xPlus_h=Get_xPlus_h(x, x_h); //Matrix of xPlus_h and values of x_h
    const Mat1 fconst=(Instnc->*f)(x, data); //Change to Mat1
    int NoOfRows=fconst.n_rows;
    int NoOfColumns=fconst.n_cols;
    Mat1 df_dx(NoOfRows, NoOfColumns*x.n_rows); //Change to Mat1
    cout<<"Jacobian NoOfRows= "<<NoOfRows<<" NoOfColumns*x.n_rows= "<<NoOfColumns*x.n_rows;
    for (int i=0; i<x.n_rows; i++)
    {
        Mat1 fx1=(Instnc->*f)(xPlus_h.row(i), data)-fconst;
        Mat1 dfx1_dx=fx1/x_h(i);
        for (int j=0; j<NoOfColumns; j++)
        {
            df_dx.col(NoOfColumns*j+i)=dfx1_dx.col(j);
        }
    }
    return df_dx;
}

inline mat Get_xPlus_h(vec x, vec& x_h)
{
    //vec x=vectorise(x1); //Ensures a column vector
    x_h=sqrt(eps(x))%x; //find x_h
    mat a;
    a<<1.0<<endr;
    double eps1=sqrt(eps(a).at(0,0));
    x_h.replace(0.0,eps1); //find elements where x_h=0 and replace them.
    //x_h.replace(0.0,.000000002); //fix if the above lines do not give good results.
    mat xPlus_h=kron(ones(x.n_rows,1),x.t())+diagmat(x_h); //Creates a matrix [x+x_h, y, z; x, y+y_h, z; x, y, z+z_h]
    return xPlus_h;
}


#endif // JACOBIAN_HPP
