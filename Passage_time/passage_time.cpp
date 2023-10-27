#include "passage_time.h"
#include <iostream>
#include <vector>
#include <cstdio>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <numeric>
passage::passage(){}
/*
 * The idea for the calculation:
 * We record each crossing events of trajectory corresponding index.
 * We need to use two vectors to do this, which is upcrossing_index and downcrossing_index;
 * To distinguish with two states, set state 0 is in the downside, and state 1 is upside.
 * Also same for the crossing events, which cross_state 0 is under local minima of the potential.
 * Here to avoid the initialize the size of vector, since we don't know how many cases will
 * happen, therefore we use two temperary vector to push them back.
 */

double cross_time(double traj_c, double traj_p, int index, double time_interval)
{
        // This function used to calculate the exact time of our passage crossing time.
        // The reason is that we record the index not exactly when crossing happends,
        // it will be always a bit more than exact time.
        double velocity = abs(traj_c-traj_p)/time_interval;
        // Calculate the velocity before pass the boundary
        double t_between_boundary = abs(abs(traj_c)-1)/velocity;
        //Calculate the absolute time
        return index*time_interval-t_between_boundary;
}


void passage::passage_time(double time_interval, double upboundary, double downboundary)
{
        //define the state for declare which states in double basin.
        double state;
        //Define how many times transfer from one state to another state.
        unsigned int jump = 0;
        int cross_state = 1;

        int beginning_state;
        int beginning_point;
        //Because we need to decide the initial state, this just to
        //initialize the vector
        for(unsigned int i; i<traj.size(); i++){
                if(traj[i] < downboundary)
                {
                        state = downboundary;
                        beginning_state = 0;
                        cross_state = 0;
                        break;
                }
                if(traj[i] > upboundary)
                {
                        state = upboundary;
                        beginning_state = 1;
                        cross_state = 1;
                        break;
                }
        }
        //For temperary saving vector, then push_back.
        std::vector<int> uptemp;
        std::cout<<"initial state is "<<beginning_state<<std::endl;
        std::vector<int> downtemp;

        for(unsigned int i=beginning_point; i<traj.size(); i++)
        {
                if(state == downboundary)
                {
                        //If the state flips, should clear the temp vector and record the first crossing
                        //events
                        if(traj[i]>=upboundary){
                                //jump to another satate
                                state = upboundary;
	                        cross_state = 1;
                                jump += 1;
                                uptemp.push_back(i);
                                // For passage time test
                                time_passage_up.push_back(i);

                                downcrossing_index.push_back(downtemp);
                                downtemp.clear();
                        }
                        else
                        {
                                if(traj[i]<downboundary & cross_state == 1)
                                {
                                        cross_state = 0;
                                        downtemp.push_back(i);
                                        time_passage_down.push_back(i);
                                }
                                if(traj[i]>downboundary & cross_state == 0 )
                                {
                                        cross_state = 1;
                                        downtemp.push_back(i);


                                        time_passage_down.push_back(i);
                                }
                        }

                }
                else
                {
                        if(traj[i]<=downboundary){
                                state = downboundary;
                                cross_state = 0;
                                jump += 1;
                                downtemp.push_back(i);

                                time_passage_down.push_back(i);

                                upcrossing_index.push_back(uptemp);
                                uptemp.clear();

                        }
                        else
                        {
                                if(traj[i]<upboundary & cross_state == 1)
                                {
                                        cross_state = 0;
                                        uptemp.push_back(i);
                                        time_passage_up.push_back(i);
                                }
                                if(traj[i]>upboundary & cross_state == 0 )
                                {
                                        cross_state = 1;
                                        uptemp.push_back(i);


                                        time_passage_up.push_back(i);
                                }
                        }
                }



        }
        std::cout<<"State index recording end!"<<std::endl;
        //For passage time calculation.
        // 1, For first-first passage time calculation
        std::cout<<"Jump time is "<<jump<<std::endl;
        std::cout<<"upcrossing_index size is "<<upcrossing_index.size()<<std::endl;

        std::cout<<"downcrossing_index size is "<<downcrossing_index.size()<<std::endl;
        
        double ffpt;
        double passage_time_temp;
        double jump_temp;
        double passage_temp;
        double transition_path_temp;
        /*
         * The idea for this part:
         * Since we got each index of the 
         * ----|-------------------|----upcrossing_index
         *     .\                  |
         *     . \     /\    /\    |
         * -------\---/--\--/--\--/------downcrossing_index
         *     .  .\ /.  .\/.  .\/..
         * corresponding index to position
         *     .  .   .  .  .  .  ..
         *    up[0]   .  .  .  .  .up[1]
         *        .   .  .  .  .  .  
         *   down[0]  1  2  3  4  5
         * The passage time in the lower well will be each index minus the next first element in
         * upcrossing_index.
         */
        for(unsigned int i=0; i<jump/2 -1; i++)
        {
                if(beginning_state == 1)
                {

                        if(upcrossing_index[i].size()>1)
                        {
                                jump_temp = cross_time(traj[downcrossing_index[i][0]], traj[downcrossing_index[i][0]-1], downcrossing_index[i][0], time_interval);
                                passage_temp = cross_time(traj[upcrossing_index[i][0]], traj[upcrossing_index[i][0]-1], upcrossing_index[i][0], time_interval);
                                ffpt = (downcrossing_index[i][0]-upcrossing_index[i][0])*time_interval;
                                // time_down =
                                FFPT.push_back(ffpt);
                                passage_time_temp = 0.;
                                transition_path_temp = jump_temp;
                                transition_path_temp-= cross_time(traj[upcrossing_index[i].back()], traj[upcrossing_index[i].back()-1], upcrossing_index[i].back(), time_interval);
                                TRAN.push_back(transition_path_temp);
                                for(size_t j=0; j<upcrossing_index[i].size()-1;j++)
                                {
                                        // passage_time_temp += (downcrossing_index[i][0]-upcrossing_index[i][j])*time_interval;
                                        passage_temp = cross_time(traj[upcrossing_index[i][j]], traj[upcrossing_index[i][j]-1], upcrossing_index[i][j], time_interval);
                                        passage_time_temp = jump_temp - passage_temp; 
                                        MFPT.push_back(passage_time_temp);
                                }
                        }
                        if(downcrossing_index[i].size()>1)
                        {
                                jump_temp = cross_time(traj[upcrossing_index[i+1][0]], traj[upcrossing_index[i+1][0]-1], upcrossing_index[i+1][0], time_interval);
                                passage_temp = cross_time(traj[downcrossing_index[i][0]], traj[downcrossing_index[i][0]-1], downcrossing_index[i][0],time_interval);
                                ffpt = jump_temp - passage_temp;
                                FFPT.push_back(ffpt);
                                passage_time_temp = 0.;
                                transition_path_temp = jump_temp;
                                transition_path_temp-= cross_time(traj[downcrossing_index[i].back()], traj[downcrossing_index[i].back()-1], downcrossing_index[i].back(), time_interval);
                                TRAN.push_back(transition_path_temp);
                                for(size_t j=0; j<downcrossing_index[i].size()-1;j++)
                                {
                                        passage_temp = cross_time(traj[downcrossing_index[i][j]], traj[downcrossing_index[i][j]-1], downcrossing_index[i][j],time_interval);
                                        passage_time_temp = jump_temp - passage_temp; 
                                        MFPT.push_back(passage_time_temp);
                                }
                        }
                }
                else
                {
                        if(downcrossing_index[i].size()>1)
                        {
                                jump_temp = cross_time(traj[upcrossing_index[i][0]], traj[upcrossing_index[i][0]-1], upcrossing_index[i][0], time_interval);
                                passage_temp = cross_time(traj[downcrossing_index[i][0]], traj[downcrossing_index[i][0]-1], downcrossing_index[i][0], time_interval);
                                ffpt = jump_temp - passage_temp;
                                FFPT.push_back(ffpt);
                                transition_path_temp = jump_temp;
                                transition_path_temp-= cross_time(traj[downcrossing_index[i].back()], traj[downcrossing_index[i].back()-1], downcrossing_index[i].back(), time_interval);
                                TRAN.push_back(transition_path_temp);
                                passage_time_temp = 0.;
                                for(size_t j=0; j<downcrossing_index[i].size();j++)
                                {
                                        passage_temp = cross_time(traj[downcrossing_index[i][j]], traj[downcrossing_index[i][j]-1], downcrossing_index[i][j], time_interval);
                                        time_passage_down.push_back(passage_temp);
                                        passage_time_temp = jump_temp - passage_temp;
                                        MFPT.push_back(passage_time_temp);
                                }
                        }
                        if(upcrossing_index[i].size()>1)
                        {
                                jump_temp = cross_time(traj[downcrossing_index[i+1][0]], traj[downcrossing_index[i+1][0]-1], downcrossing_index[i+1][0], time_interval);
                                passage_temp = cross_time(traj[upcrossing_index[i][0]], traj[upcrossing_index[i][0]-1], upcrossing_index[i][0],time_interval);
                                ffpt = jump_temp - passage_temp;
                                FFPT.push_back(ffpt);
                                transition_path_temp = jump_temp;
                                transition_path_temp -= cross_time(traj[upcrossing_index[i].back()], traj[upcrossing_index[i].back()-1], upcrossing_index[i].back(), time_interval);
                                TRAN.push_back(transition_path_temp);
                                passage_time_temp = 0.;
                                for(size_t j=0; j<upcrossing_index[i].size();j++)
                                {
                                        passage_temp = cross_time(traj[upcrossing_index[i][j]], traj[upcrossing_index[i][j]-1], upcrossing_index[i][j],time_interval);
                                        time_passage_up.push_back(passage_temp);
                                        passage_time_temp = jump_temp - passage_temp;
                                        MFPT.push_back(passage_time_temp);
                                }
                        }
                }
                if(upcrossing_index[i].size()-1>0)
                {
                	crossing_number.push_back(upcrossing_index[i].size()-1);
                	crossing_number.push_back(downcrossing_index[i].size()-1);
                	}
        }
        std::cout<<"Transfer index to time"<<std::endl;

        // For calculate the mean all first passage time;
        std::cout<<"MFPT size is "<<MFPT.size()<<std::endl;
        std::cout<<"Boundary is "<<upboundary<<" "<<downboundary;
        std::cout<<"Transition path size is "<<TRAN.size()<<std::endl;
        std::cout<<std::endl;

        std::cout<<std::endl;
        std::cout<<std::endl;

}

void passage::histogram_build(int discretize, std::vector<double> &trajectory)
{
    histogram_discretize = discretize;
    hist_plot.clear();
    hist_accum.clear();
    hist_bin.clear();
    double max_boundary = *std::max_element(trajectory.begin(), trajectory.end());
    double min_boundary = *std::min_element(trajectory.begin(), trajectory.end());
    std::vector<int> histogram_init;
    for(unsigned int i = 0; i < trajectory.size();i++ )
    {
        histogram_init.push_back(static_cast<int>(std::floor(trajectory[i]*histogram_discretize)));
    }
    int hist_min = std::floor(min_boundary*histogram_discretize);
    int hist_max = std::floor(max_boundary*histogram_discretize);
    hist_nbins = hist_max - hist_min+1;
    for(unsigned int i=0; i<hist_nbins; i++)
    {
        hist_accum.push_back(0);
        hist_bin.push_back(0);
    }
    std::cout<<"The histogram maximum:"<<hist_max/(double)histogram_discretize<<std::endl;
    hist_maximum_edge = hist_max/(double)histogram_discretize;
    std::cout<<"The number of bins is "<<hist_nbins<<std::endl;
    std::cout<<"The bins' size is "<<1/(double)discretize<<std::endl;
    //add a pseudo point in the end, max+1, the reason because we use floor function,

    for(unsigned int i = 0; i<hist_nbins; i++)
        hist_plot.push_back(0);
    for(int i = 0; i<hist_nbins; i++)
    {
        hist_plot[i]=(hist_min+i+0.5)/(double)discretize;
        if(i==0)
        {
            std::cout<<"Test: the minmum edge of histogram is "<<hist_min/(double)histogram_discretize<<std::endl;
            hist_minimum_edge = hist_min/(double)histogram_discretize;
        }
    }
    //This gives each point probability weight.
    int hist_index = 0;
    for(unsigned int i= 0; i < histogram_init.size(); i++)
    {
        hist_accum[histogram_init[i]-hist_min] +=  1;
        hist_bin[histogram_init[i]-hist_min] +=  1;
    }
    histogram_init.clear();
       //Test out:
    double sum = 0;
    for( unsigned int i = 0; i<hist_accum.size(); i++)
    {
        sum += hist_accum[i]/discretize;
    }
    // static double rescale =  sum*discretize;
    for(unsigned int i =0; i< hist_accum.size(); i++)
    {
        hist_accum[i] = hist_accum[i]/(sum*1/discretize); 
    }
    histogram_init.clear();
    sum = 0;
    for( unsigned int i = 0; i<hist_accum.size(); i++)
    {
        sum += hist_accum[i]/discretize;
    }
 


    std::cout<<"Integrate over the histogram is "<<sum<<std::endl;
    std::cout<<"size is "<<hist_plot.size()<<std::endl;
}

