/*
 * File name: gsoa.cc
 * Date:      2016/12/07 08:32
 * Author:    Jan Faigl
 */

#include <limits>
#include <chrono>
#include <boost/foreach.hpp>

#include <crl/random.h>
#include <crl/logging.h>
#include <crl/assert.h>
#include <crl/file_utils.h>
#include <crl/timerN.h>

#include "gsoa.h"
#include "gsoa_learning.h"
#include "route_path_utils.h"

using namespace std::chrono;
using namespace std;

using namespace crl;
using namespace gsoa;

#define foreach BOOST_FOREACH

typedef std::vector<int> IntVector;

/// ----------------------------------------------------------------------------
static void createPermutation(int number, IntVector &permutation)
{
   permutation.clear();
   for (int i = 0; i < number; i++)
   {
      permutation.push_back(i);
   }
}

/// ----------------------------------------------------------------------------
static void permute(IntVector &permutation)
{
   int k, tmp;
   crl::CRandom::randomize();
   for (int i = permutation.size(); i > 0; --i)
   {
      k = crl::CRandom::random() % i;
      tmp = permutation[i - 1];
      permutation[i - 1] = permutation[k];
      permutation[k] = tmp;
   }
}

/// - static method ------------------------------------------------------------
crl::CConfig &CGSOA::getConfig(crl::CConfig &config)
{
   // basic properties not included in the crl algorithm
   Base::getConfig(config);
   config.add<bool>("save-info", "disable/enable save info", true);
   config.add<bool>("save-settings", "disable/enable save settings", true);

   config.add<std::string>("result-path", "file name for the final found path (ring) as sequence of points",
                           "path");
   config.add<std::string>("result-path-dir", "dir name for the final found paths (ring) as sequence of points",
                           "paths");
   config.add<std::string>("result-tour", "file name for the final found sequence",
                           "tour");
   config.add<std::string>("result-tour-dir", "dir name for the final found sequence",
                           "tours");
   //
   // GSOA - learning parameters
   Schema::getConfig(config);
   config.add<bool>("2opt-post", "enable 2opt postprocessing of the found path", false);
   config.add<double>("radius-decrease", "Distance to decreased the radius to ensure the waypoint is in the neighbourhood", 0.1);
   config.add<bool>("in-neigh-waypoint", "enable/disable determining in neighborhood waypoint", false);
   config.add<bool>("post-opt", "enable postprocessing of the found path", false);

   //
   // Problem specification
   config.add<std::string>("problem", "Problem file");
   config.add<int>("depot-idx", "If >= 0, the particular goal is considered as the depot with the radius 0", -1);
   config.add<bool>("variable-radius", "If enabled, the input file is considered as x y radius ", false);
   config.add<double>("communication-radius", "Radius within other sensors can be read, disabled if <= 0", 0);
   config.add<std::string>("method", "Specify method in the result log", "gsoa");
   config.add<double>("sphere-degree", "Degree omege of the sphere", 10.0);

   // LKH/GLKH - parameters
   config.add<bool>("lkh", "post optimize lkh solution", false);
   config.add<bool>("glkh", "post optimize glkh solution", false);
   config.add<std::string>("glkh-instance", "path to glkh instance", "");
   config.add<int>("samples", "Number of samples in the GTSP-based instance", 1);

   return config;
}

/// - constructor --------------------------------------------------------------
CGSOA::CGSOA(crl::CConfig &config) : Base(config, "TRIAL"),
                                     DEPOT_IDX(config.get<int>("depot-idx")),
                                     VARIABLE_RADIUS(config.get<bool>("variable-radius")),
                                     SAVE_RESULTS(config.get<bool>("save-results")),
                                     SAVE_SETTINGS(config.get<bool>("save-settings")),
                                     SAVE_INFO(config.get<bool>("save-info")),
                                     IN_NEIGH_WAYPOINT(config.get<bool>("in-neigh-waypoint")),
                                     SAMPLES(config.get<int>("samples"))
{

   method = config.get<std::string>("method");

   if (!schema)
   {
      schema = new Schema(config);
   }

   const std::string fname = config.get<std::string>("problem");

   double angle, size = 0.0;
   double radius = load_goals(fname, targets, &angle, &size);

   if (name.size() == 0)
   {
      std::string n = getBasename(fname);
      size_t i = n.rfind(".tsp");
      if (i != std::string::npos)
      {
         name = n.erase(i, 4);
      }
   }

   ring = new CRing(targets, config.get<double>("radius-decrease"), IN_NEIGH_WAYPOINT);
   ring->set_center(angle, size, radius);
   ring->m = 0;
}

/// - destructor ---------------------------------------------------------------
CGSOA::~CGSOA()
{
   delete ring;
   foreach (STarget *target, targets)
   {
      delete target;
   }
}

/// - public method ------------------------------------------------------------
std::string CGSOA::getVersion(void)
{
   return "GSOA TSPNS 1.0";
}

/// - public method ------------------------------------------------------------
std::string CGSOA::getRevision(void)
{
   return "$Id: gsoa.cc 234 2018-08-16 20:09:49Z jf $";
}

/// - public method ------------------------------------------------------------
void CGSOA::visualize(const CoordsVector &path)
{
}

/// - public method ------------------------------------------------------------
void CGSOA::solve(void)
{
   crl::CRandom::randomize();
   Base::solve();
}

/// - protected method ---------------------------------------------------------
void CGSOA::load(void)
{
   DEBUG("GSOA::load -- done");
}

/// - protected method ---------------------------------------------------------
void CGSOA::initialize(void)
{
   permutation.clear();
   foreach (const STarget *target, targets)
   {
      permutation.push_back(target->label);
   }
}

/// - protected method ---------------------------------------------------------
void CGSOA::after_init(void)
{
}

/// - protected method ---------------------------------------------------------
void CGSOA::iterate(int iter)
{
   crl::CRandom::randomize();

   finalSolution.clear();

   int finalBestSolutionStep;

   TargetPtrVector allTargets;

   CoordsVector pts;
   for (auto target : targets)
   {
      pts.push_back(target->coords);
   }
   ring->initialize_neurons(pts);

   foreach (STarget *target, targets)
   {
      target->selectedWinner = 0;
      target->stepWinnerSelected = -1;
   }
   schema->G = config.get<double>("learning-gain");
   schema->mi = config.get<double>("learning-rate");
   //thresholds for the termination conditions
   const double MAX_ERROR = config.get<double>("termination-error");
   const int MAX_STEPS = config.get<int>("termination-max-steps");
   const bool TERM_CHANGE = config.get<bool>("termination-change");

   double error = 2 * MAX_ERROR;
   int step = 0;
   CoordsVector solution;
   CoordsVector bestSolution;
   IntVector sequence;
   int bestSolutionStep = -1;
   double bestSolutionLength = std::numeric_limits<double>::max();
   const bool BEST_SOLUTION = config.get<bool>("best-solution");
   bool term = false;
   IntVector routes[2];
   int routeCur = 0;
   int routePrev = 1;
   while (!((error < MAX_ERROR)) && (step < MAX_STEPS) && not term)
   { //perform adaptation step
      error = refine(step, error);

      getSolution(step, solution); //collect solution

      const double length = get_path_length(solution, ring->radius);

      if (length < bestSolutionLength)
      {
         bestSolution = solution;
         bestSolutionStep = step;
         bestSolutionLength = length;
      }

      schema->G = schema->G * (1 - schema->GAIN_DECREASING_RATE * (step + 1));

      step++;
      
   } //end step loop

   tSolve.stop();

   fillResultRecord(iter);

   resultLog
       << bestSolutionLength //
       << bestSolutionLength //
       << step
       << step
       << finalBestSolutionStep
       << config.get<bool>("post-opt")
       << config.get<bool>("2opt-post")
       << bestSolutionLength
       << 0
       << MAX_STEPS
       << schema->GAIN_DECREASING_RATE
       << config.get<double>("learning-gain")
       << crl::result::endrec;
   DEBUG("Best solution with the length: " << bestSolutionLength << " found in: " << bestSolutionStep << " steps");
}

/// - protected method ---------------------------------------------------------
double CGSOA::refine(int step, double errorMax)
{
   double errorToGoal = errorMax;
   double error = 0.0;
   permute(permutation);

   schema->updateExp(targets.size(), step);
   for (IntVector::iterator i = permutation.begin(); i != permutation.end(); i++)
   {
      STarget *target = targets[*i];
      SNeuron *prevWinner = target->stepWinnerSelected == step - 1 ? target->selectedWinner : 0;
      int r = 0;
      SWinnerSelection *winner = ring->selectWinner(step, target, errorToGoal);
      if (winner and winner->hasWinner)
      {
         ring->adapt(step);
         if (error < errorToGoal)
         {
            error = errorToGoal; //update error
         }
      }
   } //end permutation of all targets
   ring->regenerate(step);

   return error; // return largest error to city
}

/// - protected method ---------------------------------------------------------
void CGSOA::save(void)
{
   std::string dir;
   updateResultRecordTimes(); //update timers as load and initilization is outside class
   if (SAVE_SETTINGS)
   {
      saveSettings(getOutputIterPath(config.get<std::string>("settings"), dir));
   }
   if (SAVE_INFO)
   {
      saveInfo(getOutputIterPath(config.get<std::string>("info"), dir));
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
void CGSOA::release(void)
{
}

/// - protected method ---------------------------------------------------------
void CGSOA::defineResultLog(void)
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
      resultLog << result::newcol << "BEST_LENGTH";
      resultLog << result::newcol << "STEPS";
      resultLog << result::newcol << "CURRENT_STEP";
      resultLog << result::newcol << "SOLUTION_STEP";
      resultLog << result::newcol << "POST_OPT";
      resultLog << result::newcol << "2OPT_POST";
      resultLog << result::newcol << "LENGTH";
      resultLog << result::newcol << "SAMPLES";
      resultLog << result::newcol << "MAX_STEPS";
      resultLog << result::newcol << "GAIN_RATE";
      resultLog << result::newcol << "GAIN";
      resultLogInitialized = true;
   }
}

/// - protected method ---------------------------------------------------------
void CGSOA::fillResultRecord(int trial)
{
   resultLog << result::newrec << name << method << trial;
   long t[3] = {0l, 0l, 0l};
   tLoad.addTime(t);
   tInit.addTime(t);
   tSolve.addTime(t);
   tSave.addTime(t);
   std::string tmp = std::to_string((float)lkhtime);
   resultLog << t[0] << t[1] << t[2] << tmp << (opttime);
}

/// - private method -----------------------------------------------------------
void CGSOA::drawPath(void)
{
}

/// - private method -----------------------------------------------------------
void CGSOA::drawRing(int step)
{
}

/// - private method -----------------------------------------------------------
void CGSOA::savePic(int step, bool detail, const std::string &dir_suffix)
{
   static int lastStep = step;
   static int i = 0;
   if (lastStep != step)
   {
      i = 0;
   }
   if (canvas)
   {
      canvas->redraw();
      std::string dir;
      std::string file = getOutputIterPath(config.get<std::string>("pic-dir") + dir_suffix, dir);
      assert_io(createDirectory(file), "Cannot create file in path '" + file + "'");
      std::stringstream ss;
      ss << file << "/"
         << "iter-" << std::setw(3) << std::setfill('0') << step;
      ss << "-" << std::setw(4) << std::setfill('0') << i;

      std::string suffixes(config.get<std::string>("pic-ext"));
      if (!suffixes.empty())
      {
         std::string::size_type cur = 0;
         std::string::size_type next;
         do
         {
            next = suffixes.find(',', cur);
            const std::string &ext = suffixes.substr(cur, next - cur);
            if (!ext.empty())
            {
               assert_io(canvas->save(ss.str() + "." + ext), "Cannot create output canvas file '" + file + "'");
            }
            cur = next + 1;
         } while (next != std::string::npos);
      }
      else
      {
         ss << "." << config.get<std::string>("pic-ext");
         assert_io(canvas->save(ss.str()), "Cannot create output canvas file '" + ss.str() + "'");
      }
   }
   lastStep = step;
   i++;
}

/// - private method -----------------------------------------------------------
void CGSOA::getSolution(int step, CoordsVector &solution) const
{
   ring->get_ring_path(step, solution);
}

/// - private method -----------------------------------------------------------
bool CGSOA::testSolution(void) const
{
   // bool ok = true;
   int ok = 0;

   for (auto &target : targets)
   {
      for (auto &pt : finalSolution)
      {
         if (pt.angle(target->coords) <= ring->omega + 0.0001)
         {
            ok += 1;
            break;
         }
      }
   }
   return ok == targets.size();
}
/// - private method -----------------------------------------------------------
void CGSOA::saveSolution(std::string prefix) const
{
}

/// - private method -----------------------------------------------------------
void CGSOA::postOptimizeSolution(CoordsVector &route)
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

      Coords act2 = closest_with_optimization(ring->omega, prev, next, target->coords, tmp);
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
/* end of gsoa.cc */
