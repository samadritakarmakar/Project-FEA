#ifndef CHECKSUPPORTEDORDER_HPP
#define CHECKSUPPORTEDORDER_HPP
#include <iostream>
void CheckSupportedOrder(const int order, const int orderLimit)
{
    if (order>orderLimit)
    {
        std::cout<<"Given Order is Greater than the Order Limit of "<<orderLimit<<"!\n";
        throw;
    }
}


#endif // CHECKSUPPORTEDORDER_HPP
