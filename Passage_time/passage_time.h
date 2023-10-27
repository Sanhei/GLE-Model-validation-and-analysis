#ifndef PASSAGE
#define PASSAGE

#include <vector>
using namespace std;

class passage
{
public:
        passage();
        double time_interval;
        vector<double> traj;

        vector<vector<int>> upcrossing_index;
        vector<vector<int>> downcrossing_index;

        vector<double> MFPT;
        vector<double> FFPT;
        vector<double> TRAN;
        vector<vector<double>> all_passage_time;
        vector<int> crossing_number;



        void trajectory_input(vector<double> traj);
        void passage_time(double time_interval, double upboundary, double downboundary);
        void histogram_build(int discretize, std::vector<double> &trajectory);

        // For histogram
        vector<double> hist_plot;
        vector<double> hist_accum;
        vector<double> hist_bin;

        vector<double> time_passage_up;
        vector<double> time_passage_down;




        double hist_minimum_edge;
        double hist_maximum_edge;
        int histogram_discretize;
        int hist_nbins;

        double mfpt;
        double ffpt;

};
#endif
