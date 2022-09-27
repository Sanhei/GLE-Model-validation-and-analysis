#ifndef MEMORY
#define MEMORY

#include <vector>
#include <math.h>

class memory
{
public:
    memory();
    std::vector<double> G;
    std::vector<double> Gsmooth;
    std::vector<double> kernel;
    void memory_G_t(double time_interval, std::vector<double> corvv, std::vector<double> corux, int memory_range);
    void memory_kernel();
private:
    double time_interval_;
};
#endif
