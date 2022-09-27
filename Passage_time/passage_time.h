#ifndef PASSAGE
#define PASSAGE

#include <vector>
using namespace std;

class passage
{
public:
        passage();
        double upboundary;
        double downboundary;
        double time_interval;
        vector<vector<int>> upcrossing_index;
        vector<vector<int>> downcrossing_index;

        vector<double> MFPT;
        vector<double> FFPT;
        vector<vector<double>> all_passage_time;
        vector<int> crossing_number;


        void passage_time(vector<double> traj, double time_interval);
        void histogram_build(int discretize, std::vector<double> &trajectory);
        vector<double> hist_plot;
        vector<double> hist_accum;
        double hist_minimum_edge;
        double hist_maximum_edge;
        int histogram_discretize;
        int hist_nbins;
};
#endif
