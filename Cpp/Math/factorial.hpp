#ifndef FACTORIAL_HPP
#define FACTORIAL_HPP
#include <armadillo>

unsigned int factorial(unsigned int n) 
{ 
    int res = 1, i; 
    for (i = 2; i <= n; i++) 
        res *= i; 
    return res; 
} 



#endif
