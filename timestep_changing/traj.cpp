#include <iostream>
#include "parameter.h"
//#include "potential.hpp"
#include "v_verlet.hpp"
#include <fstream>
#include <vector>
#include <string>
using namespace std;
int main()
{
    double x = -1.0; // initial position
    double temp=0;
    double xp = x;
    double y  = 0;
    //std::vector<double> tlist = {0.01, 0.001, 0.0001};
    
    //std::vector<double> tlist = {0.0075};
    std::cout<<"tau_D is "<<tau_D<<std::endl;
    std::cout<<"tau is "<<tau<<std::endl;
    std::cout<<"tau_m is "<<tau_m<<std::endl;

    double h=0.001;
    //to generate` different files
    string base(".txt");

    std::cout<<"Simulation begins"<<std::endl;

    std::cout<<"For time interval "<<h<<" simulation starts"<<std::endl;
    ofstream file;
    std::string filename(".txt");
    file.open(std::to_string(h)+filename);
    //file << "# time   position"<< std::endl;
    //file << "#" << std::endl;
    for(size_t i=0; i<length; i++){
        temp = x;
        x = x_update(x, xp, y, h);
        y = y_update(y, temp, h);
        xp = temp;
	file << x << std::endl;
    }
    file.close();
}
 
