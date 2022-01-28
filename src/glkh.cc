/*
 * File name: glkh.cc
 * Date:      2020/12/07 08:32
 * Author:    Jindriska Deckerova
 */

#include <limits>
#include <chrono>
#include <boost/foreach.hpp>

#include <crl/random.h>
#include <crl/logging.h>
#include <crl/assert.h>
#include <crl/file_utils.h>
#include <crl/timerN.h>

#include "glkh.h"
#include "route_path_utils.h"

using namespace chrono;
using namespace std;

using namespace crl;
using namespace glkh;

#define foreach BOOST_FOREACH

/// - static method ------------------------------------------------------------
static void process_first_comment(const string &line, double *angle, double *size)
{
    //COMMENT: type = sphere, neigh_ang = 0.0872664625997, neigh_size = 10
    string::size_type sz;
    string delimiter = ":", delimiter1 = "=";
    size_t pos = line.find(delimiter);
    string data = line.substr(0 + line.substr(0, pos).length(), line.length());

    vector<string> sublines;
    string token;
    delimiter = ",";
    while ((pos = data.find(delimiter)) != string::npos)
    {
        token = data.substr(0, pos);
        sublines.push_back(string(token));
        data.erase(0, pos + delimiter.length());
    }
    sublines.push_back(data);

    string angl = sublines[1].substr(sublines[1].find(delimiter1) + delimiter1.length(), sublines[1].length());
    *angle = stod(angl, &sz);

    string siz = sublines[2].substr(sublines[2].find(delimiter1) + delimiter1.length(), sublines[2].length());
    *size = stod(siz, &sz);
}

/// - static method ------------------------------------------------------------
static void process_second_comment(const string &line, double *radius)
{
    //COMMENT: sphere_limit = 1.0471975512, rect_limit = [-1.0471975512, 0.785398163397, -0.785398163397, 1.57079632679]

    string::size_type sz;
    string delimiter = ":", delimiter1 = "=";
    size_t pos = line.find(delimiter);
    string data = line.substr(0 + line.substr(0, pos).length(), line.length());

    vector<string> sublines;
    string token;
    delimiter = ",";
    while ((pos = data.find(delimiter)) != string::npos)
    {
        token = data.substr(0, pos);
        sublines.push_back(string(token));
        data.erase(0, pos + delimiter.length());
    }
    sublines.push_back(data);

    string rad = sublines[0].substr(sublines[0].find(delimiter1) + delimiter1.length(), sublines[0].length());
    *radius = stod(rad, &sz);
}

/// - static method ------------------------------------------------------------
crl::CConfig &CGLKH::getConfig(crl::CConfig &config)
{
    return config;
}

/// - constructor --------------------------------------------------------------
CGLKH::CGLKH(crl::CConfig &config) : Base(config, "TRIAL"),
                                     DEPOT_IDX(config.get<int>("depot-idx")),
                                     VARIABLE_RADIUS(config.get<bool>("variable-radius")),
                                     SAVE_RESULTS(config.get<bool>("save-results")),
                                     SAVE_SETTINGS(config.get<bool>("save-settings")),
                                     SAVE_INFO(config.get<bool>("save-info")),
                                     SAMPLES(config.get<int>("samples"))
{

    method = config.get<string>("method");
    const string fname = config.get<string>("problem");

    radius = load(fname, targets, &omega, &neigh_size);
    samples = config.get<int>("samples");

    if (name.size() == 0)
    {
        string n = getBasename(fname);
        size_t i = n.rfind(".tsp");
        if (i != string::npos)
        {
            name = n.erase(i, 4);
        }
    }
}

/// - destructor ---------------------------------------------------------------
CGLKH::~CGLKH()
{
    foreach (gsoa::STarget *target, targets)
    {
        delete target;
    }
}

/// - public method ------------------------------------------------------------
string CGLKH::getVersion(void)
{
    return "GLKH TSPNS 1.0";
}

/// - public method ------------------------------------------------------------
string CGLKH::getRevision(void)
{
    return "$Id: glkh.cc 234 2018-08-16 20:09:49Z jf $";
}

/// - public method ------------------------------------------------------------
void CGLKH::solve(void)
{
    crl::CRandom::randomize();
    Base::solve();
}

/// - protected method ---------------------------------------------------------
double CGLKH::load(const string &filename, gsoa::TargetPtrVector &targets, double *angle, double *size)
{
    Coords pt;
    ifstream in(filename.c_str());
    string line;
    stringstream ss;
    bool loading = false, neig_info = false, sphere_info = false;
    int idx = 0;
    double lat, lon;
    double radius = 0.0;
    *angle = *size = 0.0;

    // COMMENT: type = band, neigh_ang = 0.0872664625997, neigh_size = 10
    // COMMENT: sphere_limit = 1.0471975512, rect_limit = [-1.0471975512, 0.785398163397, -0.785398163397, 1.57079632679]
    string delimiter = ":";
    string delimiter1 = "=";
    string::size_type sz;
    while (getline(in, line))
    {
        if (line == "NODE_COORD_SECTION:")
        {
            loading = true;
            continue;
        }
        if (line == "GTSP_SET_SECTION:")
        {
            break;
        }
        string comment = line.substr(0, line.find(delimiter));
        if (comment == "COMMENT")
        {
            if (!neig_info)
            {
                neig_info = true;
                process_first_comment(line, angle, size);
            }
            else if (!sphere_info)
            {
                sphere_info = true;
                process_second_comment(line, &radius);
            }
        }

        if (loading)
        {
            ss.clear();
            ss << line;
            ss >> idx >> lon >> lat;
            lat = lat * (M_PI / 180);
            lon = lon * (M_PI / 180);

            radius = 1.0;
            double x = radius * cos(lon) * cos(lat); //cos(lat) * sin(lon);
            double y = radius * sin(lon) * cos(lat); //sin(lat) * sin(lon);
            double z = radius * sin(lat);            //cos(lon);
            double radius = radius;

            targets.push_back(new gsoa::STarget(targets.size(), Coords(x, y, z, radius), *angle, *size));
        }
    }

    DEBUG("GLKH::load -- done");
    return radius;
}

/// - protected method ---------------------------------------------------------
void CGLKH::initialize(void)
{
    if (config.get<bool>("lkh"))
    {
        method = "lkh";
    }
    else if (config.get<bool>("glkh"))
    {
        method = "glkh";
    }
}

/// - protected method ---------------------------------------------------------
void CGLKH::after_init(void)
{
}

/// - protected method ---------------------------------------------------------
void CGLKH::load(void)
{
}

/// - protected method ---------------------------------------------------------
void CGLKH::iterate(int iter)
{
    string problem;
    if (config.get<bool>("glkh"))
    {
        problem = config.get<string>("glkh-instance") ;
    } else {
        problem = config.get<string>("problem");
    }
    compute_tour(problem, iter);

    crl::CRandom::randomize();    

    finalSolution.clear();
    finalSequence.clear();

    if (config.get<bool>("lkh"))
    {
        cout << method << endl;
        method = "lkh";
        load_lkh("glkh/GLKH-1.0/tmp.tour", targets, finalSolution, finalSequence);
    }
    else if (config.get<bool>("glkh"))
    {
        method = "glkh";
        load_glkh("glkh/GLKH-1.0/tmp.tour", problem,
                  targets, finalSolution, finalSequence, samples, radius, omega);
    }

    double length = get_path_length(finalSolution, radius);
    double optlength = length;
    if (config.get<bool>("post-opt"))
    {
        auto start = high_resolution_clock::now();

        postOptimizeSolution(finalSolution);
        optlength = get_path_length(finalSolution, radius);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);

        stringstream ss;
        ss << duration.count();
        opttime = stod(ss.str());
        DEBUG("Length: " << length << " after post-opt: " << optlength);
    }

    fillResultRecord(iter);

    resultLog
        << optlength //
        << -1
        << -1
        << config.get<bool>("post-opt")
        << 0 // config.get<bool>("2opt-post")
        << length
        << samples
        << crl::result::endrec;
    DEBUG("Best solution with the length: " << length);
}

/// - protected method ---------------------------------------------------------
void CGLKH::compute_tour(const string &problem, const int iter)
{
    system("rm glkh/GLKH-1.0/tmp.par");
    system("rm glkh/GLKH-1.0/tmp.tsp");
    system("rm glkh/GLKH-1.0/tmp.tour");

    stringstream ss;
    ss << "cp " << problem << " glkh/GLKH-1.0/tmp.tsp";
    cout << ss.str() << endl;
    system(ss.str().c_str());

    int seed = (int) crl::CRandom::random(1, 100);
    string text = "";
    text += "PROBLEM_FILE = tmp.tsp\n";
    text += "RUNS = 10\n";                      // TODO:
    text += "SEED = " + to_string(seed) + "\n"; // TODO:
    text += "OUTPUT_TOUR_FILE = tmp.tour\n";

    ofstream ofs("glkh/GLKH-1.0/tmp.par");
    ofs << text;
    ofs.close();

    // Run GLKH solver
    auto start = high_resolution_clock::now();
    system("cd glkh/GLKH-1.0 && ./GLKH tmp.par");
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    ss.str("");
    ss << duration.count();
    lkhtime = stod(ss.str());

}

/// - protected method ---------------------------------------------------------
void CGLKH::save(void)
{
    string dir;
    updateResultRecordTimes(); //update timers as load and initilization is outside class
    if (SAVE_SETTINGS)
    {
        saveSettings(getOutputIterPath(config.get<string>("settings"), dir));
    }
    if (SAVE_INFO)
    {
        saveInfo(getOutputIterPath(config.get<string>("info"), dir));
    }
    if (SAVE_RESULTS)
    {
        string file = getOutputIterPath(config.get<string>("result-path"), dir);
        assert_io(createDirectory(dir), "Can not create file in the path'" + file + "'");

        const int i = 0;
        stringstream ss;
        ss << file << "-" << setw(2) << setfill('0') << i << ".txt";
        ofstream ofs(ss.str());
        assert_io(not ofs.fail(), "Cannot create path '" + ss.str() + "'");
        ofs << setprecision(14);
        foreach (const Coords &pt, finalSolution)
        {
            ofs << pt.x << " " << pt.y << " " << pt.z << endl;
        }
        assert_io(not ofs.fail(), "Error occur during path saving");
        ofs.close();

        file = getOutputIterPath(config.get<string>("result-tour"), dir);
        assert_io(createDirectory(dir), "Can not create file in the path'" + file + "'");

        ss.str("");
        ss << file << "-" << setw(2) << setfill('0') << i << ".txt";
        ofs.open(ss.str());
        assert_io(not ofs.fail(), "Cannot create path '" + ss.str() + "'");
        ofs << setprecision(14);
        foreach (const int s, finalSequence)
        {
            ofs << s << endl;
        }
        assert_io(not ofs.fail(), "Error occur during path saving");
        ofs.close();
    }
}

/// - protected method ---------------------------------------------------------
void CGLKH::release(void)
{
}

/// - protected method ---------------------------------------------------------
void CGLKH::defineResultLog(void)
{
    static bool resultLogInitialized = false;
    if (!resultLogInitialized)
    {
        resultLog << result::newcol << "NAME";
        resultLog << result::newcol << "METHOD";
        resultLog << result::newcol << "TRIAL";
        resultLog << result::newcol << "RTIME";
        resultLog << result::newcol << "CTIME";
        resultLog << result::newcol << "UTIME";
        resultLog << result::newcol << "LKHTIME";
        resultLog << result::newcol << "OPTTIME";
        resultLog << result::newcol << "OPT_LENGTH";
        resultLog << result::newcol << "STEPS";
        resultLog << result::newcol << "SOLUTION_STEP";
        resultLog << result::newcol << "POST_OPT";
        resultLog << result::newcol << "2OPT_POST";
        resultLog << result::newcol << "LENGTH";
        resultLog << result::newcol << "SAMPLES";
        resultLogInitialized = true;
    }
}

/// - protected method ---------------------------------------------------------
void CGLKH::fillResultRecord(int trial)
{
    resultLog << result::newrec << name << method << trial;
    long t[3] = {0l, 0l, 0l};
    tLoad.addTime(t);
    tInit.addTime(t);
    tSolve.addTime(t);
    tSave.addTime(t);
    string tmp = to_string((float)lkhtime);
    resultLog << t[0] << t[1] << t[2] << tmp << (opttime);
}

/// - private method -----------------------------------------------------------
void CGLKH::getSolution(int step, CoordsVector &solution) const
{
    // solution.clear();
    // for (auto target : finalSolution)
    // {
    //     double x = target.x;
    //     double y = target.y;
    //     double z = target.z;
    //     double lat = atan2(z, (sqrt(x * x + y * y))) * (180 / M_PI);
    //     double lon = atan2(y, x) * (180 / M_PI);
    //     solution.push_back(Coords(lat, lon));
    // }
}

/// - private method -----------------------------------------------------------
bool CGLKH::testSolution(void) const
{
    // bool ok = true;
    int ok = 0;

    for (auto &target : targets)
    {
        for (auto &pt : finalSolution)
        {
            if (pt.angle(target->coords) <= omega + 0.0001)
            {
                ok += 1;
                break;
            }
        }
    }
    return ok == targets.size();
}
/// - private method -----------------------------------------------------------
void CGLKH::saveSolution(string prefix) const
{
    const double final_length = get_path_length(finalSolution, radius, true);
    string fname = "solutions/" + prefix + "_" + name + "_" + to_string(final_length) + "_.sol";
    ofstream ofs(fname);

    for (auto target : finalSolution)
    {
        double x = target.x;
        double y = target.y;
        double z = target.z;
        double lat = atan2(z, (sqrt(x * x + y * y))) * (180 / M_PI);
        double lon = atan2(y, x) * (180 / M_PI);
        ofs << target.x << " " << target.y << " " << target.z << "\n";
    }
    ofs.close();
}

/// - private method -----------------------------------------------------------
void CGLKH::postOptimizeSolution(CoordsVector &route)
{
    Coords x;
    int n = finalSequence.size();
    for (int i = 0; i < 5 * n; i++)
    {
        auto &target = targets[finalSequence[i % n]];

        auto &prev = route[(i - 1 + n) % n];
        auto &act = route[(i + 0) % n];
        auto &next = route[(i + 1) % n];
        Coords tmp;

        Coords act2 = closest_with_optimization(omega, prev, next, target->coords, tmp);
        if (prev.angle(act2) + next.angle(act2) <= prev.angle(act) + next.angle(act))
        {
            act = act2;
        }
        else
        {
            if (prev.angle(act2) + next.angle(act2) > prev.angle(act) + next.angle(act) + 1e-5)
            {
                ERROR("This should not happen -> improve closest 2 by local optimization");
            }
        }
    }
}

/* end of glkh.cc */
