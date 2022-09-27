/*
 * Copyright © 2022-2023 Qingyuan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include "memory_extraction/potential.h"
#include "memory_extraction/correlation.h"
#include <cmath>
#include "memory_extraction/spline.h"
#include <fftw3.h>
#include "figure_plot.h"
#include <typeinfo>
#include "memory_extraction/memory.h"
//For spline
#include <boost/math/interpolators/makima.hpp>

//For automatically get time_interval from filename.
#include <experimental/filesystem>
#include <iterator>
#include <algorithm>
#include "Passage_time/passage_time.h"

void analytical_potential(std::vector<double> &traj, std::vector<double> &potential_a)
{
    double depth = 1.0;
    for(unsigned int i; i<traj.size(); i++){
        potential_a.push_back(depth*(4*pow(traj[i], 3) - 4*traj[i]));
    }
}
namespace fs = std::experimental::filesystem;
int main(int argc, char **argv)
{
    double depth = 1.0;
    //This step to get the time interval, since we name our file by time steps.
        // The start path. Use Program path
        const fs::path startPath{ fs::path("./timestep_changing") };

        // Here we will store all file names
        std::vector<fs::path> files{};
        std::string ss;
        std::copy_if(fs::directory_iterator(startPath), {}, std::back_inserter(files), [](const fs::directory_entry& de) { return de.path().extension() == ".txt";});
        std::cout << files[0].string().substr(20, 21) << '\n';
        ss = files[0].string().substr(20, 21);
        std::cout<<ss.substr(0, 8);






    double time_interval = std::stod(ss.substr(0, 8))*10;
    std::cout<<"Time interval is "<<time_interval<<std::endl;
    double temperature   = 1.0;
    std::ifstream inputfile;
    inputfile.open(files[0]);
    if(!inputfile)
    {
        std::cout<<"unable to openfile";
        exit(1);
    }
    double x, y;
    //Since in the trajectory one column is time, another is position;
    std::vector<double> trajectory;
    std::cout<<"Reading trajectory"<<std::endl;
    while(inputfile>>x)
    {
        trajectory.push_back(x);
    }





	std::cout<<"Calculating the number of crossing events"<<std::endl;


        //Passage calculation
        passage Passage_cnum_calc;
        Passage_cnum_calc.passage_time(trajectory, time_interval);
        int discretize = 1;
	std::vector<double> number_of_crossing(Passage_cnum_calc.crossing_number.begin(), Passage_cnum_calc.crossing_number.end());
        Passage_cnum_calc.histogram_build(discretize, number_of_crossing);
        number_of_crossing.clear();
        Passage_cnum_calc.crossing_number.clear();
        // 1. Record the histogram of number_of_crossing;
        std::ofstream hist_num_cross;
        hist_num_cross.open("./figure/hist_num_cross.txt");
        for(size_t i=0; i<Passage_cnum_calc.hist_plot.size(); i++)
        {
                hist_num_cross<<Passage_cnum_calc.hist_plot[i]<<" "<<Passage_cnum_calc.hist_accum[i]<<std::endl;
        }
        hist_num_cross.close();
        //Clear the vector for next histogram calculation;
        Passage_cnum_calc.hist_plot.clear();
        Passage_cnum_calc.hist_accum.clear();

        std::cout<<"Calculating Mean first first passage time"<<std::endl;


        // 2. For Mean first-first passage time calculation
        //We need to change this for Mean first passage time calculation
        discretize = 5;
        Passage_cnum_calc.histogram_build(discretize, Passage_cnum_calc.MFPT);
        std::ofstream hist_MFPT;
        hist_MFPT.open("./figure/hist_MFPT.txt");
        for(size_t i=0; i<Passage_cnum_calc.hist_plot.size(); i++)
        {
                hist_MFPT<<Passage_cnum_calc.hist_plot[i]<<" "<<Passage_cnum_calc.hist_accum[i]<<std::endl;
        }
        hist_MFPT.close();

        std::cout<<"Calculating Mean first second passage time"<<std::endl;
        // Release memory
        //
        Passage_cnum_calc.hist_plot.clear();
        Passage_cnum_calc.hist_accum.clear();
        Passage_cnum_calc.MFPT.clear();


        Passage_cnum_calc.histogram_build(discretize, Passage_cnum_calc.FFPT);
        std::ofstream hist_FFPT;
        hist_FFPT.open("./figure/hist_FFPT.txt");
        for(size_t i=0; i<Passage_cnum_calc.hist_plot.size(); i++)
        {
                hist_FFPT<<Passage_cnum_calc.hist_plot[i]<<" "<<Passage_cnum_calc.hist_accum[i]<<std::endl;
        }
        hist_FFPT.close();

        Passage_cnum_calc.hist_plot.clear();
        Passage_cnum_calc.hist_accum.clear();
        Passage_cnum_calc.FFPT.clear();



    int distcretize = 50;
    histogram_constructor Hist_Obj;
    Hist_Obj.SetTemperature(temperature);
    Hist_Obj.histogram_build(distcretize, trajectory);
    Hist_Obj.boundary_influence(-1.5, 1.6);
    //Here we set the boundary is -1.5 and 1.5;
    std::vector<double> Free_energy;
    for(unsigned int i; i<Hist_Obj.hist_accum.size(); i++)
    {
        Free_energy.push_back(-temperature*log(Hist_Obj.hist_accum[i]));
    }
    std::cout<<"Histogram end"<<std::endl;
    std::cout<<"The boundary influnce"<<std::endl;
    /*
     * In this part we will integrate over the boundary.
     * This case, our well position is 1 and -1, so we set
     * the boundary is -1.5 and 1.5.
     * The method to study is to integrate the area beyond
     * this part([-1.5, 1.5]) over the whole integrate.
     * So we can see the rate of this part in the trajectory.
     */





     //For theoretical plot
    std::vector<double> theoretical;
    std::cout<<"Size is "<<Hist_Obj.hist_plot.size()<<std::endl;
    std::vector<double> plot_edge;

    double shift; //for error analyse, make the simulation result and theoretical overlap
    for(unsigned int i=0; i < Hist_Obj.hist_plot.size(); i++)
    {
        plot_edge.push_back(Hist_Obj.hist_plot[i]);
        theoretical.push_back(depth*(pow(Hist_Obj.hist_plot[i],4)-2*pow(Hist_Obj.hist_plot[i], 2))+3);
        if(i == Hist_Obj.midpoint_index)
                shift = theoretical[i]-Hist_Obj.free_energy[i];
    }
    // Do the shift and see the result
    Hist_Obj.shift_free_energy(shift);



//Get the function of Free energy( spline )
    /* Two spline mthods:
     * 1. Cubic spline, use 3 order polynomiers ax^3+bx^2+cx+d;
     *    tk::spline, need to call spline.h file.
     * 2. Akima spline, the local suddenly drop or raise doesn't
     *    change the other distinct, which means, the function has
     *    less vibration;
     */
    tk::spline Potential;
    Potential.set_points(Hist_Obj.hist_plot, Free_energy, tk::spline::cspline);
    //using boost::math::interpolators::makima;
    //auto Potential = makima<std::vector<double>>(std::move(Hist_Obj.hist_plot), std::move(Hist_Obj.free_energy));
    //Spline function saved in Potential;
       //For spline plot
    int splinesize = int((Hist_Obj.hist_maximum_edge - Hist_Obj.hist_minimum_edge)/(double)0.0001);
    //std::cout<<"Spline fault"<<std::endl;
    std::vector<double> smallbins;
    std::vector<double> Energy;
    for(unsigned int i=0; i < splinesize; i++)
    {
        smallbins.push_back(Hist_Obj.hist_minimum_edge+0.0001*i);
        Energy.push_back(Potential(Hist_Obj.hist_minimum_edge+0.0001*i));
    }
    /* This will plot the potential, and smooth function
     * smallbins is the for spline function, show the detail of
     * spline result.
     */
    //Figure.figure_plot::potential_plot(plot_edge,smallbins, Energy, Free_energy, theoretical, Hist_Obj.free_energy_shift);

    //Clear the memory of shift free energy;
    Hist_Obj.free_energy_shift.clear();


    //std::cout<<typeid(Potential).name()<<" is the type of Potential"<<std::endl;



    //Gradient of potentail
    std::cout<<"Gradient of potentail"<<std::endl;
    std::vector<double> gradient_u_plot;
    std::vector<double> analytical_gradient;
    std::vector<double> gradient_spline_plot;
    Hist_Obj.derivative(plot_edge, Free_energy, gradient_u_plot);

    analytical_potential(plot_edge, analytical_gradient);


    for(unsigned int i=0; i<plot_edge.size(); i++){
        gradient_spline_plot.push_back(plot_edge[i]);
    }
    //Spline;
    //using boost::math::interpolators::makima;
    //auto G_u = makima<std::vector<double>>(std::move(plot_edge), std::move(gradient_u_plot));
    tk::spline G_u;
    G_u.set_points(plot_edge, gradient_u_plot, tk::spline::cspline);

    std::cout<<"spline gradient right"<<std::endl;
    std::vector<double> gradient_u_spline;
    for(unsigned int i=0; i<gradient_spline_plot.size(); i++){
        if(gradient_spline_plot[i]>2)
                std::cout<<"This number is rong"<<std::endl;
        gradient_u_spline.push_back(G_u(gradient_spline_plot[i]));
    }
    std::cout<<"right"<<std::endl;
    //Figure.figure_plot::gradient_plot(gradient_spline_plot, gradient_u_spline, analytical_gradient);


    //Plot gradien
    //Release memory
    Free_energy.clear();
    theoretical.clear();
    plot_edge.clear();
    smallbins.clear();
    Energy.clear();
    gradient_u_plot.clear();
    analytical_gradient.clear();
    gradient_spline_plot.clear();
    gradient_u_spline.clear();

    Hist_Obj.gradient_plot.clear();
    Hist_Obj.gradient_U.clear();
    Hist_Obj.free_energy_shift.clear();










    //Correlation Calculation
    std::cout<<"Correlation of velocity_and velocity"<<std::endl;
    correlation_calculation Cor_vv;
    Cor_vv.N = trajectory.size() - 2;
    Cor_vv.SetTime_interval(time_interval);
    Cor_vv.correlation_calculation::velocity(trajectory);
    Cor_vv.correlation_calculation::fft(Cor_vv.velocity_, Cor_vv.velocity_);
    Cor_vv.correlation_calculation::Convolution(Cor_vv.velocity_, Cor_vv.velocity_);
    //Figure.figure_plot::Corvv_plot(time_interval, Cor_vv.Correlation);
    std::ofstream corvv_r;
    corvv_r.open("./figure/corvv.txt");
    for(size_t i; i<Cor_vv.Correlation.size();++i)
            corvv_r<<Cor_vv.Correlation[i]<<std::endl;
    corvv_r.close();

    Cor_vv.Correlation.clear();







    //End of the  Corvv and begin the Cor_gux
    std::cout<<"The correlation of gradient Potential and position"<<std::endl;

    correlation_calculation Cor_gux;
    Cor_gux.SetTime_interval(time_interval);
    Cor_gux.N = trajectory.size() - 2;
    //Transfer the vector to fftw_complex
    static fftw_complex *trajectory_f;
    trajectory_f = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(trajectory.size()-2));
    for(unsigned int i = 0; i<trajectory.size()-2; i++)
    {
        trajectory_f[i][0] = trajectory[i+1];
    }
    static fftw_complex *gradient_u;
    gradient_u = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(trajectory.size()-2));
    for(unsigned int i=0; i<trajectory.size()-2; i++)
    {
        gradient_u[i][0] = G_u(trajectory[i+1]);
        //gradient_u[i][0] = depth*(4*pow(trajectory[i], 3) - 4*trajectory[i]);
        //gradient_u[i][0] = spline(trajectory[i]);
    }
    trajectory.clear();
    Cor_gux.correlation_calculation::fft(trajectory_f, trajectory_f);
    Cor_gux.correlation_calculation::fft(gradient_u, gradient_u);
    Cor_gux.correlation_calculation::Convolution(gradient_u, trajectory_f);
    //Figure.figure_plot::Corg_ux_plot(time_interval, Cor_gux.Correlation);
    std::ofstream corgux_r;
    std::cout<<"Recording Corux"<<std::endl;
    corgux_r.open("./figure/corux.txt");
    for(size_t i; i<Cor_gux.Correlation.size();++i)
            corgux_r<<Cor_gux.Correlation[i]<<std::endl;
    corgux_r.close();
    //Figure.figure_plot::Correlation_plot(time_interval, Cor_gux.Correlation, Cor_vv.Correlation);
    std::cout<<"Untill"<<std::endl;


    std::ifstream input_corvv;
    input_corvv.open("./figure/corvv.txt");
    if(!input_corvv)
    {
        std::cout<<"Unable to open Corvv"<<std::endl;
        exit(1);
    }
    std::vector<double> corvv;
    std::cout<<"Reading Corvv"<<std::endl;
    while(input_corvv>>x)
    {
        corvv.push_back(x);
    }
    input_corvv.close();
    std::ifstream input_corgux;
    input_corgux.open("./figure/corux.txt");
    std::vector<double> corgux;
    if(!input_corgux)
    {
        std::cout<<"Unable to open Corgux"<<std::endl;
        exit(1);
    }
    std::cout<<"Reading corgux"<<std::endl;
    while(input_corgux>>x)
    {
        corgux.push_back(x);
    }
    input_corgux.close();
    //Figure.Corg_ux_plot(time_interval, corgux);
    std::cout<<"Subplot"<<std::endl;
    //This is plot for two correlation, the end 1000000 is the plot range, since
    //Figure.Correlation_plot(time_interval, corvv, corgux, 100000);
    //our memory not enough for plotting them both.
    std::cout<<"Memory calculating"<<std::endl;
    memory Memory;
    Memory.memory_G_t(time_interval, corvv, corgux, 2000);
    std::cout<<"Image plotting"<<std::endl;
    std::cout<<"happend"<<std::endl;

    //Our G(t) plotting behaviour is like vibration, here we do the average;
    Memory.memory_kernel();
    corvv.clear();
    corgux.clear();
    std::ofstream memory_output;
    memory_output.open("kernel.txt");
    for(unsigned int i = 0; i<Memory.kernel.size(); i++)
            memory_output<<time_interval*i*2<< "  "<<Memory.kernel[i]<< std::endl;
    //Figure.memory_G_t_kernel_plot(time_interval, Memory.Gsmooth, Memory.kernel);


    std::string sr = "fitting.py";
    std::string command = "python ";
    command += sr;
    system(command.c_str());

    std::cout<<"END"<<std::endl;
    return 0;

}
