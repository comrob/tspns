/*
 * File name: glkh.h
 * Date:      2020/12/07 08:32
 * Author:    Jindriska Deckerova
 */

#ifndef __GLKH_H__
#define __GLKH_H__

#include <crl/config.h>
#include <crl/alg/algorithm.h>

#include "coords.h"
#include "target.h"

using namespace std;

typedef vector<int> IntVector;

static void load_lkh(const string &filename, const gsoa::TargetPtrVector &targets, CoordsVector &solution, IntVector &sequence)
{
    ifstream in(filename.c_str());
    string line;
    stringstream ss;
    bool loading = false;
    int idx = 0, seq = 0;

    while (getline(in, line))
    {
        if (line == "-1")
        {
            break;
        }
        if (idx == 6)
        {
            loading = true;
        }
        if (loading)
        {
            ss.clear();
            ss << line;
            ss >> seq;
            solution.push_back(targets[seq - 1]->coords);
            sequence.push_back(seq - 1);
        }
        idx++;
    }
}
static void load_glkh(const string &filename, const string &gtspfilename, const gsoa::TargetPtrVector &targets, CoordsVector &solution, IntVector &sequence, const int &m, const double &radius, const double &omega)
{
    ifstream in(filename.c_str());
    string line;
    stringstream ss;
    bool loading = false;
    int idx = 0, seq = 0, sam = 0, sample;
    CoordsVector samples;

    for (int i = 0; i < m; i++)
    {
        samples.push_back(Coords(radius * cos(omega), radius * sin(omega) * cos((i * 2 * M_PI) / m), radius * sin(omega) * sin((i * 2 * M_PI) / m)).norm2());
    }

    for (auto &target : targets)
    {
        target->samples.clear();
        double theta = -atan2(target->coords.z, sqrt(target->coords.x * target->coords.x + target->coords.y * target->coords.y));
        double phi = atan2(target->coords.y, target->coords.x);
        for (int i = 0; i < m; i++)
        { // rotate using rotation matrix with shuffled rows (2., 3., 1.) 
            auto pt = samples[i];
            target->samples.push_back(Coords(pt.x * cos(theta) * cos(phi) + (-1) * pt.y * sin(phi) + pt.z * sin(theta) * cos(phi),
                                             pt.x * cos(theta) * sin(phi) + pt.y * cos(phi) + pt.z * sin(theta) * sin(phi),
                                             (-1) * pt.x * sin(theta) + 0 + pt.z * cos(theta)));
            target->samples[i] = target->samples[i].norm2();
        }
    }

    while (getline(in, line))
    {
        if (line == "-1")
        {
            break;
        }
        if (idx == 6)
        {
            loading = true;
        }
        if (loading)
        {
            ss.clear();
            ss << line;
            ss >> sample;
            seq = (sample / m >= targets.size() || sample % m == 0) ? (sample / m) - 1 : sample / m;
            sam = (sample % m) > 0 ? (sample % m) - 1 : m - 1;
            solution.push_back(targets[seq]->samples[sam]);
            sequence.push_back(seq);
        }
        idx++;
    }
}

namespace glkh
{

    class CGLKH : public crl::CAlgorithm
    {
        typedef crl::CAlgorithm Base;

    public:
        static crl::CConfig &getConfig(crl::CConfig &config);

        CGLKH(crl::CConfig &config);
        ~CGLKH();

        string getVersion(void);
        string getRevision(void);

        void solve(void);

    protected:
        void load();
        double load(const string &filename, gsoa::TargetPtrVector &targets, double *angle, double *size);
        void initialize(void);
        void after_init(void);

        void iterate(int iter);
        void compute_tour(const string &problem, const int iter);

        void save(void);
        void release(void);

        void defineResultLog(void);
        void fillResultRecord(int trial);

    private:
        void getSolution(int step, CoordsVector &solution) const;
        bool testSolution(void) const;
        void saveSolution(string prefix = "glkh") const;
        void postOptimizeSolution(CoordsVector &route);

    private:
        const int DEPOT_IDX;
        const bool VARIABLE_RADIUS;

        const bool SAVE_RESULTS;
        const bool SAVE_SETTINGS;
        const bool SAVE_INFO;

        const int SAMPLES;

        string method;
        double lkhtime;
        double opttime;

        gsoa::TargetPtrVector targets;
        CoordsVector finalSolution;
        IntVector finalSequence;

        int samples;
        double radius;
        double omega;
        double neigh_size;
    };

} // end name glkh

#endif

/* end of glkh.h */
