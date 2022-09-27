#ifndef FIGURE_PLOT
#define FIGURE_PLOT
#include <vector>

class figure_plot
{
public:
    figure_plot();
    void potential_plot(std::vector<double> hist_plot, std::vector<double> smallbins, std::vector<double> Spline, std::vector<double> free_energy, std::vector<double> theoretical, std::vector<double> shift_free_energy);
    void gradient_plot(std::vector<double> x, std::vector<double> gradient_spline, std::vector<double> y);
    void Corvv_plot(double timestep, std::vector<double> Corvv);
    void Corg_ux_plot(double timestep, std::vector<double> Corg_ux);
    void Correlation_plot(double timestep, std::vector<double> Corvv, std::vector<double> Corg_ux, int plot_range);
    void single_plot(std::vector<double> x, std::vector<double> y);
    void memory_G_t_kernel_plot(double timestep, std::vector<double> G_t, std::vector<double> kernel);
};
#endif
