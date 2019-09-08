#ifndef LOCALINTEGRATION_HPP
#define LOCALINTEGRATION_HPP
#include "Form.hpp"
#include "FormMultiThread.hpp"
#include "TrialFunction.hpp"
#include "TrialFunctionNeumannSurface.hpp"
#include "TrialFunctionNeumannLine.hpp"
#include "TestFunctionGalerkin.hpp"
template <class GenericTrialFunction>
class LocalIntegrator
{
public:
    LocalIntegrator(Form<GenericTrialFunction>& a, GenericTrialFunction& u,
                    TestFunctionGalerkin<GenericTrialFunction>& v):
        u(u),v(v)
    {
        _a.SetNumOfThreads(1);
        _a.form[0]=a;
        ResultingMat=std::vector<sp_mat>(1);
        ResultingVector=std::vector<mat>(1);
        ResultingScalar=std::vector<double>(1);
        usingFormMultiThread = false;
    }

    LocalIntegrator(FormMultiThread<GenericTrialFunction>& a, GenericTrialFunction& u,
                    TestFunctionGalerkin<GenericTrialFunction>& v):
        _a(a),u(u),v(v)
    {
        numOfThreads= _a.GetNumOfThreads();
        ResultingMat=std::vector<sp_mat>(numOfThreads);
        ResultingVector=std::vector<mat>(numOfThreads);
        ResultingScalar=std::vector<double>(numOfThreads);
        usingFormMultiThread = true;
    }

    virtual ~LocalIntegrator()
    {
        //----Do-Nothing
    }

    virtual sp_mat weak_form(Form<GenericTrialFunction>& a, GenericTrialFunction& u,
                             TestFunctionGalerkin<GenericTrialFunction>& v)
    {
        cout<<"Warining!!! Internal default virtual function for matrix in a single thread is running\n";
        return a.inner(a.grad(v),a.grad(u))*a.dX(u);
    }

    virtual sp_mat weak_form(FormMultiThread<GenericTrialFunction>& a, GenericTrialFunction& u,
                     TestFunctionGalerkin<GenericTrialFunction>& v, int thread)
    {
        cout<<"Warining!!! Internal default virtual function for matrix in multi thread is running\n";
        return a[thread].inner(a[thread].grad(v),a[thread].grad(u))*a[thread].dX(u);
    }

    virtual mat weak_form_vector(Form<GenericTrialFunction>& a, GenericTrialFunction& u,
                                 TestFunctionGalerkin<GenericTrialFunction>& v)
    {
        vec b;
        cout<<"Warining!!! Internal default virtual function for vector in a single thread is running\n";
        b<<0<<endr<<0<<endr<<-9.81<<endr;
        return a.dot(v,b)*a.dX(u);
    }

    virtual mat weak_form_vector(FormMultiThread<GenericTrialFunction>& a, GenericTrialFunction& u,
                                 TestFunctionGalerkin<GenericTrialFunction>& v, int thread)
    {
        vec b;
        b<<0<<endr<<0<<endr<<-9.81<<endr;
        cout<<"Warining!!! Internal default virtual function for matrix in multi thread is running\n";
        return a[thread].dot(v,b)*a[thread].dX(u);
    }

    virtual double scalar_integration(Form<GenericTrialFunction>& a, GenericTrialFunction& u)
    {
        return a.dX(u);
    }

    virtual double scalar_integration(FormMultiThread<GenericTrialFunction>& a, GenericTrialFunction& u, int thread)
    {
        return a[thread].dX(u);
    }


    void local_intergrator()
    {
        _a[0].set_u_Internal(u);
        _a[0].set_v_Internal(v);
        //cout<<"Element Number ="<<_a[0].ElementNumber<<" Gauss pt is at "<<_a[0].GaussPntr<<"\n";
        ResultingMat[0]=weak_form(_a[0],u,v);
        //cout<<mat(a.ResultingMat)<<"\n";
        int NoOfGaussPnts= GetNumOfGaussPoints(0);
        //cout<<"No of GaussPnts = "<<NoOfGaussPnts<<"\n";
        for (int i=0; i<NoOfGaussPnts-1; i++)
        {
            _a[0].NextGaussPntr();
            //cout<<"Element Number ="<<_a[0].ElementNumber<<"Gauss pt is at "<<_a[0].GaussPntr<<" i= "<<i<<"\n";
            ResultingMat[0]=ResultingMat[0]+weak_form(_a[0],u,v);
            //if(_a[0].GaussPntr==215)
                //cout<<"Resulting Mat size = "<<_a[0].ResultingMat.n_rows<<", "<<_a[0].ResultingMat.n_cols<<"\n";
        }
        _a[0].GaussPntr=0;
    }

    void local_intergrator(int thread)
    {
        _a[thread].set_u_Internal(u);
        _a[thread].set_v_Internal(v);
        //cout<<"Element Number ="<<_a[0].ElementNumber<<" Gauss pt is at "<<_a[0].GaussPntr<<"\n";
        ResultingMat[thread]=weak_form(_a,u,v, thread);
        //cout<<mat(a.ResultingMat)<<"\n";
        int NoOfGaussPnts=GetNumOfGaussPoints(thread);
        //cout<<"No of GaussPnts = "<<NoOfGaussPnts<<"\n";
        for (int i=0; i<NoOfGaussPnts-1; i++)
        {
            _a[thread].NextGaussPntr();
            //cout<<"Element Number ="<<_a[0].ElementNumber<<"Gauss pt is at "<<_a[0].GaussPntr<<" i= "<<i<<"\n";
            ResultingMat[thread]=ResultingMat[thread]+weak_form(_a,u,v, thread);
            //if(_a[0].GaussPntr==215)
                //cout<<"Resulting Mat size = "<<_a[0].ResultingMat.n_rows<<", "<<_a[0].ResultingMat.n_cols<<"\n";
        }
        _a[thread].GaussPntr=0;
    }

    void local_intergrator_vector()
    {
        _a[0].set_u_Internal(u);
        _a[0].set_v_Internal(v);
        //cout<<"Element Type ="<<a.ElementType<<"Gauss pt is at "<<a.GaussPntr<<"\n";
        ResultingVector[0]=weak_form_vector(_a[0],u,v);
        int NoOfGaussPnts=GetNumOfGaussPoints(0);
        for (int i=0; i<NoOfGaussPnts-1; i++)
        {
            _a[0].NextGaussPntr();
            //cout<<"Element Type ="<<a.ElementType<<"Gauss pt is at "<<a.GaussPntr<<"i= "<<i<<"\n";
            ResultingVector[0]=ResultingVector[0]+weak_form_vector(_a[0],u,v);
        }
        _a[0].GaussPntr=0;
    }

    void local_intergrator_vector(int thread)
    {
        _a[thread].set_u_Internal(u);
        _a[thread].set_v_Internal(v);
        //cout<<"Element Type ="<<a.ElementType<<"Gauss pt is at "<<a.GaussPntr<<"\n";
        ResultingVector[thread]=weak_form_vector(_a,u,v, thread);
        int NoOfGaussPnts=GetNumOfGaussPoints(thread);
        for (int i=0; i<NoOfGaussPnts-1; i++)
        {
            _a[thread].NextGaussPntr();
            //cout<<"Element Type ="<<a.ElementType<<"Gauss pt is at "<<a.GaussPntr<<"i= "<<i<<"\n";
            ResultingVector[thread]=ResultingVector[thread]+weak_form_vector(_a,u,v, thread);
        }
        _a[thread].GaussPntr=0;
    }

    void local_integration_scalar()
    {
        _a[0].set_u_Internal(u);
        _a[0].set_v_Internal(v);
        //cout<<"Element Type ="<<_a[0].ElementType<<"Gauss pt is at "<<_a[0].GaussPntr<<"\n";
        ResultingScalar[0]=scalar_integration(_a[0],u);
        int NoOfGaussPnts=GetNumOfGaussPoints(0);
        for (int i=0; i<NoOfGaussPnts-1; i++)
        {
            _a[0].NextGaussPntr();
            //cout<<"Element Type ="<<_a[0].ElementType<<"Gauss pt is at "<<_a[0].GaussPntr<<"\n";
            ResultingScalar[0]=ResultingScalar[0]+scalar_integration(_a[0],u);
        }
        _a[0].GaussPntr=0;
    }

    void local_integration_scalar(int thread)
    {
        _a[thread].set_u_Internal(u);
        _a[thread].set_v_Internal(v);
        //cout<<"Element Type ="<<_a[thread].ElementType<<"Gauss pt is at "<<_a[thread].GaussPntr<<"\n";
        ResultingScalar[thread]=scalar_integration(_a,u, thread);
        int NoOfGaussPnts=GetNumOfGaussPoints(thread);
        for (int i=0; i<NoOfGaussPnts-1; i++)
        {
            _a[thread].NextGaussPntr();
            //cout<<"Element Type ="<<_a[thread].ElementType<<"Gauss pt is at "<<_a[thread].GaussPntr<<"i= "<<i<<"\n";
            ResultingScalar[thread]=ResultingScalar[thread]+scalar_integration(_a,u, thread);
        }
        _a[thread].GaussPntr=0;
    }


    sp_mat GetResultingMat(int thread=0)
    {
        return ResultingMat[thread];
    }

    mat GetResultingVector(int thread=0)
    {
        return ResultingVector[thread];
    }

    double GetResultingScalar(int thread=0)
    {
        return ResultingScalar[thread];
    }

    int GetElementNumber(int thread=0)
    {
        return _a[thread].ElementNumber;
    }

    void GoToNextElement(int thread=0)
    {
        _a[thread].NextElementNumber();
    }

    void GoToNextElementType(int thread=0)
    {
        _a[thread].NextElementType();
    }

    bool IsUsingFormMultiThread()
    {
        return usingFormMultiThread;
    }

    void SetElementStartTo(int thread, int ElementStart)
    {
        _a[thread].SetElementStartTo(ElementStart);
    }

    int GetNumOfGaussPoints(int thread=0)
    {
        int v_numOfGaussPoints=v.numOfGaussPoints[_a[thread].ElementType];
        int u_numOfGaussPoints=u.GetNumberOfGaussPoints(_a[thread].ElementType);
        int NumOfGaussPoints=std::min(v_numOfGaussPoints, u_numOfGaussPoints);
        //cout<<"NumOfGaussPoints= "<<NumOfGaussPoints<<"\n";
        return NumOfGaussPoints;
    }

private:
    //Form<GenericTrialFunction>& a;
    FormMultiThread<GenericTrialFunction> _a;
    GenericTrialFunction& u;
    TestFunctionGalerkin<GenericTrialFunction>& v;
    std::vector<sp_mat> ResultingMat;
    std::vector<mat> ResultingVector;
    std::vector<double> ResultingScalar;
    int numOfThreads;
    bool usingFormMultiThread;
};
#endif // LOCALINTEGRATION_HPP
