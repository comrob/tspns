/*
 * File name: gsoa.h
 * Date:      2016/12/07 08:32
 * Author:    Jan Faigl
 */

#ifndef __GSOA_H__
#define __GSOA_H__

#include <crl/config.h>
#include <crl/alg/algorithm.h>

#include <crl/gui/shape.h>

#include "coords.h"
#include "target.h"
#include "gsoa_ring.h"


static void process_first_comment(const std::string &line, double *angle, double *size)
{
   //COMMENT: type = sphere, neigh_ang = 0.0872664625997, neigh_size = 10
   std::string::size_type sz;
   std::string delimiter = ":", delimiter1 = "=";
   size_t pos = line.find(delimiter);
   std::string data = line.substr(0 + line.substr(0, pos).length(), line.length());

   std::vector<std::string> sublines;
   std::string token;
   delimiter = ",";
   while ((pos = data.find(delimiter)) != std::string::npos)
   {
      token = data.substr(0, pos);
      sublines.push_back(std::string(token));
      data.erase(0, pos + delimiter.length());
   }
   sublines.push_back(data);

   std::string angl = sublines[1].substr(sublines[1].find(delimiter1) + delimiter1.length(), sublines[1].length());
   *angle = std::stod(angl, &sz);

   std::string siz = sublines[2].substr(sublines[2].find(delimiter1) + delimiter1.length(), sublines[2].length());
   *size = std::stod(siz, &sz);
}

static void process_second_comment(const std::string &line, double *radius)
{
   //COMMENT: sphere_limit = 1.0471975512, rect_limit = [-1.0471975512, 0.785398163397, -0.785398163397, 1.57079632679]

   std::string::size_type sz;
   std::string delimiter = ":", delimiter1 = "=";
   size_t pos = line.find(delimiter);
   std::string data = line.substr(0 + line.substr(0, pos).length(), line.length());

   std::vector<std::string> sublines;
   std::string token;
   delimiter = ",";
   while ((pos = data.find(delimiter)) != std::string::npos)
   {
      token = data.substr(0, pos);
      sublines.push_back(std::string(token));
      data.erase(0, pos + delimiter.length());
   }
   sublines.push_back(data);

   std::string rad = sublines[0].substr(sublines[0].find(delimiter1) + delimiter1.length(), sublines[0].length());
   *radius = std::stod(rad, &sz);
}

/// ----------------------------------------------------------------------------
static double load_goals(const std::string &filename, gsoa::TargetPtrVector &targets, double *angle, double *size)
{
   Coords pt;
   std::ifstream in(filename.c_str());
   std::string line;
   std::stringstream ss;
   bool loading = false, neig_info = false, sphere_info = false;
   int idx = 0;
   double lat, lon;
   double radius = 0.0;
   *angle = *size = 0.0;

   // COMMENT: type = band, neigh_ang = 0.0872664625997, neigh_size = 10
   // COMMENT: sphere_limit = 1.0471975512, rect_limit = [-1.0471975512, 0.785398163397, -0.785398163397, 1.57079632679]
   std::string delimiter = ":";
   std::string delimiter1 = "=";
   std::string::size_type sz;
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
      // std::cout << line << "\n";
      std::string comment = line.substr(0, line.find(delimiter));
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

   return radius;
}

// static void load_lkh(const std::string &filename, const gsoa::TargetPtrVector &targets, CoordsVector &solution, gsoa::IntVector &sequence)
// {
//    std::ifstream in(filename.c_str());
//    std::string line;
//    std::stringstream ss;
//    bool loading = false;
//    int idx = 0, seq = 0;

//    while (getline(in, line))
//    {
//       if (line == "-1")
//       {
//          break;
//       }
//       if (idx == 6)
//       {
//          loading = true;
//       }
//       if (loading)
//       {
//          ss.clear();
//          ss << line;
//          ss >> seq;
//          solution.push_back(targets[seq - 1]->coords);
//          sequence.push_back(seq - 1);
//       }
//       idx++;
//    }
// }
// static void load_glkh(const std::string &filename, const std::string &gtspfilename, const gsoa::TargetPtrVector &targets, CoordsVector &solution, gsoa::IntVector &sequence, const int &m, const double &radius, const double &omega)
// {
//    std::ifstream in(filename.c_str());
//    std::string line;
//    std::stringstream ss;
//    bool loading = false;
//    int idx = 0, seq = 0, sam = 0, sample;
//    CoordsVector samples;

//    for (int i = 0; i < m; i++)
//    {
//       samples.push_back(Coords(radius * cos(omega), radius * sin(omega) * cos((i * 2 * M_PI) / m), radius * sin(omega) * sin((i * 2 * M_PI) / m)).norm2());
//    }

//    for (auto &target : targets)
//    {
//       target->samples.clear();
//       double theta = -atan2(target->coords.z, sqrt(target->coords.x * target->coords.x + target->coords.y * target->coords.y));
//       double phi = atan2(target->coords.y, target->coords.x);
//       for (int i = 0; i < m; i++)
//       { // rotate using rotation matrix with shuffled rows (2., 3., 1.) - TODO: why shuffled?
//          auto pt = samples[i];
//          target->samples.push_back(Coords(pt.x * cos(theta) * cos(phi) + (-1) * pt.y * sin(phi) + pt.z * sin(theta) * cos(phi),
//                                           pt.x * cos(theta) * sin(phi) + pt.y * cos(phi) + pt.z * sin(theta) * sin(phi),
//                                           (-1) * pt.x * sin(theta) + 0 + pt.z * cos(theta)));
//          target->samples[i] = target->samples[i].norm2();
//       }
//    }

//    while (getline(in, line))
//    {
//       if (line == "-1")
//       {
//          break;
//       }
//       if (idx == 6)
//       {
//          loading = true;
//       }
//       if (loading)
//       {
//          ss.clear();
//          ss << line;
//          ss >> sample;
//          seq = (sample / m >= targets.size() || sample % m == 0) ? (sample/m)-1 : sample/m;
//          sam = (sample % m) > 0 ? (sample % m) - 1 : m - 1;
//          solution.push_back(targets[seq]->samples[sam]);
//          sequence.push_back(seq);
//       }
//       idx++;
//    }
// }

namespace gsoa {

   class CGSOA : public crl::CAlgorithm {
      typedef crl::CAlgorithm Base;
      typedef std::vector<int> IntVector;
      public:

      static crl::CConfig &getConfig(crl::CConfig &config);

      CGSOA(crl::CConfig &config);
      ~CGSOA();

      std::string getVersion(void);
      std::string getRevision(void);

      void visualize(const CoordsVector & path);
      void solve(void);

      protected:
      void load(void);
      void initialize(void);
      void after_init(void); 

      void iterate(int iter);
      double refine(int step, double errorMax);

      void save(void);
      void release(void);

      void defineResultLog(void);
      void fillResultRecord(int trial);

      private:
      void drawPath(void);
      void drawRing(int step);
      void savePic(int step, bool detail = false, const std::string &dir_suffix = "");

      void getSolution(int step, CoordsVector &solution) const;
      bool testSolution( void ) const;
      void saveSolution( std::string prefix = "gsoa" ) const;
      void postOptimizeSolution( CoordsVector & route );

      private:
      const int DEPOT_IDX; 
      const bool VARIABLE_RADIUS;
      const bool SAVE_RESULTS;
      const bool SAVE_SETTINGS;
      const bool SAVE_INFO;

      const bool IN_NEIGH_WAYPOINT;
      const int SAMPLES;

      std::string method;
      double lkhtime;
      double opttime;

      IntVector permutation;
      TargetPtrVector targets;
      CoordsVector finalSolution;
      IntVector finalSequence;

      CRing *ring;
   };

} // end name gsoa

#endif

/* end of gsoa.h */
