#ifndef FORMMULTITHREAD_HPP
#define FORMMULTITHREAD_HPP
#include "Form.hpp"
#include <thread>

template<class GenericTrialFunction>
class FormMultiThread
{
public:
    FormMultiThread()
    {
        std::thread detectCores;
        int numOfCores=detectCores.hardware_concurrency();
        SetNumOfThreads(numOfCores);
    }

    FormMultiThread(int NumOfThreads)
    {
        if(NumOfThreads<=0)
        {
            FormMultiThread();
        }
        else
        {
            SetNumOfThreads(NumOfThreads);
        }
    }

    void SetNumOfThreads(int NumOfThreads)
    {
        numOfThreads=NumOfThreads;
        form=std::vector<Form<GenericTrialFunction>>(numOfThreads);
    }

    Form<GenericTrialFunction>& operator[](unsigned int ThreadNum)
    {
        return form[ThreadNum];
    }

    int GetNumOfThreads()
    {
        return numOfThreads;
    }


    std::vector<Form<GenericTrialFunction>> form;
private:
    int numOfThreads;
};
#endif // FORMMULTITHREAD_HPP
