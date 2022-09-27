#include "figure_plot.h"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
figure_plot::figure_plot(){}
void figure_plot::potential_plot(std::vector<double> hist_plot,std::vector<double> smallbins,  std::vector<double> spline, std::vector<double> free_energy, std::vector<double> theoretical, std::vector<double> shift_free_energy)
{
    plt::figure();
    plt::plot( hist_plot, free_energy, {{"label", "Histogram result"}});
    std::cout<<"Histogram plot"<<std::endl;
    plt::plot( smallbins, spline, {{"label", "Spline result"}});
    std::cout<<"Spline plot"<<std::endl;
    plt::plot( hist_plot, theoretical, {{"label", "Theoretical"}});
    plt::plot( hist_plot, shift_free_energy, {{"label", "Shift of free energy"}});
    std::cout<<"Theoretical"<<std::endl;
    plt::xlabel("Position");
    plt::ylabel("Potentail");
    plt::legend();
    plt::show();
}
void figure_plot::Corvv_plot(double timestep, std::vector<double> Corvv)
{
    std::vector<double> time;
    for(int i; i<Corvv.size(); i++)
    {
        time.push_back(i*timestep);
    }
    plt::figure();
    plt::xlabel("Time");
    plt::ylabel("Correlation of velocity and velocity");
    plt::semilogx(time, Corvv);
    plt::title("Correlation");
    plt::show();
}
void figure_plot::Corg_ux_plot(double timestep, std::vector<double> Corvv)
{
    std::vector<double> time;
    for(int i; i<Corvv.size(); i++)
    {
        time.push_back(i*timestep);
    }
    plt::figure();
    plt::xlabel("Time");
    plt::ylabel("Correlation of gradient of Potentail and position");
    plt::semilogx(time, Corvv);
    plt::title("Correlation");
    plt::show();
}
template <typename T>
std::vector<T> slicing(std::vector<T> const& v, int X)
//This only for Correlation plot, since we cannot plot two figure in one graph.(it takes all the memory)
//So we cut it start from 0, and end with X.
{
        //Begin and End iterator
        auto first = v.begin();
        auto last = v.begin() + X+1;
        //Copy the elements

        std::vector<T> vector(first, last);

        //return to results;
        return vector;
}
void figure_plot::Correlation_plot(double timestep, std::vector<double> Corvv, std::vector<double> Corg_ux, int plot_range)
{
    std::vector<double> time;
    std::vector<double> corvv;
    std::vector<double> corgux;
    corvv  = slicing(Corvv,   plot_range);
    corgux = slicing(Corg_ux, plot_range);
    //The above code for select amount of elements to show since our memory not enough
    for(int i; i<corvv.size(); i++)
    {
        time.push_back(i*timestep);
    }//to plot the whole dataset.
    //Py_Initialize();
    //Calling python
    //PyRun_SimpleString("import matplotlib.pyplot as plt");
    //PyRun_SimpleString("plt.subplot(2, 2, 1)");
    //plt::suptitle("Correlation of velocity and velocity");
    const long nrows = 1, ncols = 2;
    long row = 0, col = 0;
    long spanr = 1, spanc = 2;
    plt::subplot2grid(nrows, ncols, row, col, spanr, spanr);
    plt::semilogx(time, corvv);
    plt::xlabel("time($\\tau$)");
    plt::ylabel("$C^{vv}$");
    //plt::suptitle("Correlation of gradient of potential and position");
    //PyRun_SimpleString("plt.subplot(2, 2, 2)");
    long col2 = 1;
    plt::subplot2grid(nrows, ncols, row, col2, spanr, spanr);
    plt::semilogx(time, corgux);
    //PyRun_SimpleString("plt.xlabel(r'time($\\tau$)')");
    //PyRun_SimpleString("plt.ylabel(r'$C^{\\triangledown U q}$')");
    plt::xlabel("time($\\tau$)");
    plt::ylabel("$C^{\\triangledown U x}$");
    const char* filename = "./figure/correlation.svg";
    std::cout<<"Saving correlation in figure."<<std::endl;
    plt::save(filename);
    plt::show();
    corvv.clear();
    corgux.clear();
    time.clear();
    std::cout<<"End of Correlation plotting"<<std::endl;
    //Py_Finalize();
    //Python end
}


void figure_plot::gradient_plot(std::vector<double> x, std::vector<double> gradient_spline, std::vector<double> y)
{
    plt::figure();
    plt::xlabel("Postion");
    plt::ylabel("Gradinet of potential");
    plt::plot(x, gradient_spline, {{"label", "Spline result"}});
    plt::plot(x, y, {{"label","Derivative of Potential"}});
    std::cout<<"analytical is right"<<std::endl;
    plt::legend();
    plt::show();
}
void figure_plot::single_plot(std::vector<double> x, std::vector<double> y)
{
    //Calling python
    Py_Initialize();
    std::cout<<"testing"<<std::endl;
    std::cout<<"testing"<<std::endl;
    plt::plot(x, y);
    std::cout<<"after"<<std::endl;
    plt::show();
    //Python end
    Py_Finalize();
}



void figure_plot::memory_G_t_kernel_plot(double timestep, std::vector<double> G_t, std::vector<double> kernel)
{
    std::vector<double> time;
    std::vector<double> time2;
    //The above code for select amount of elements to show since our memory not enough
    for(int i=0; i<G_t.size(); i++)
    {
        time.push_back(i*2*timestep);
    }//to plot the whole dataset.
    for(int i=0; i<kernel.size(); ++i)
        time2.push_back(i*2*timestep);
    //Py_Initialize();
    //Calling python
    //PyRun_SimpleString("import matplotlib.pyplot as plt");
    //PyRun_SimpleString("plt.subplot(2, 2, 1)");
    //plt::suptitle("Correlation of velocity and velocity");
    const long nrows = 1, ncols = 2;
    long row = 0, col = 0;
    long spanr = 1, spanc = 2;
    plt::subplot2grid(nrows, ncols, row, col, spanr, spanr);
    plt::plot(time, G_t);
    plt::xlabel("time($\\tau$)");
    plt::ylabel("$G(t)$");
    //plt::suptitle("Correlation of gradient of potential and position");
    //PyRun_SimpleString("plt.subplot(2, 2, 2)");
    long col2 = 1;
    plt::subplot2grid(nrows, ncols, row, col2, spanr, spanr);
    plt::plot(time2, kernel);
    //PyRun_SimpleString("plt.xlabel(r'time($\\tau$)')");
    //PyRun_SimpleString("plt.ylabel(r'$C^{\\triangledown U q}$')");
    plt::xlabel("time($\\tau$)");
    plt::ylabel("$Memory kernel$");
    const char* filename = "./figure/G_t_memory.svg";
    std::cout<<"Saving G(t) and Memory kernel in figure."<<std::endl;
    plt::save(filename);
 
    plt::show();
    time.clear();
    std::cout<<"End of Memory plotting"<<std::endl;
    //Py_Finalize();
    //Python end
}


