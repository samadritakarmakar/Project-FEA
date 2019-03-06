#include "UserFunctions.h"
void setConstants(Models &Model1)
{
    for (int i=0;i<Model1.NoOfConstants;i++)
    {
        std::cout<<"Enter the value of constant "<<(i+1)<<" ";
        std::cin>>Model1.constants(i);
    }
} 
