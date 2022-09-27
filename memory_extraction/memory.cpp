#include "memory.h"
#include <vector>
#include <iostream>

memory::memory(){}
void memory::memory_G_t(double time_interval, std::vector<double> corvv, std::vector<double> corux, int memory_range)
{
        time_interval_ = time_interval;
        //Pass the value for kernel calculating
    std::cout<<"Memory kernel calculating" << std::endl;
    G.push_back(0);
    G.push_back(2/(double)(time_interval*corvv[0])*(corux[1+1]-corux[0+1]/(double)corvv[0]*corvv[1]));
    double tempsum = 0;
    for(unsigned int i=2; i< memory_range;++i)
    {
        tempsum = 0.0;
        for(unsigned int j=1; j<i-1; ++j)
        {
                tempsum += G[i-j]*corvv[j];
        }
        G.push_back(2/(double)(time_interval*corvv[0])*(corux[i+1]-corux[0+1]/corvv[0]*corvv[i] - time_interval*tempsum));
    }
    std::cout<<"Memory lenth is"<<G.size();
}
void memory::memory_kernel()
{
    for(unsigned int i=0; i<G.size()/2; ++i)
    {
        Gsmooth.push_back((G[2*i+1]+G[2*i])/2);
    }
    //Calculating memory kernel gamma.
    for(unsigned int i=0; i<Gsmooth.size() - 1; ++i)
    {
        kernel.push_back((Gsmooth[i+1]-Gsmooth[i])/time_interval_/2);
    }
    G.clear();
}
