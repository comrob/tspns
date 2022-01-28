/*
 * File name: route_path_utils.cc
 * Date:      2016/12/10 18:13
 * Author:    Jan Faigl
 */

#include <cmath>
#include <cstdio>

#include "route_path_utils.h"

/// - function -----------------------------------------------------------------
double get_path_length(const CoordsVector &pts, const double R, bool closed)
{
   double len = 0.;
   for (int i = 1; i < pts.size(); ++i)
   {
   double l = (pts[i - 1].geom_distance(pts[i], R));
  // std::cout << pts[i-1] << "->" << pts[i] << " = " << l << std::endl;
      len += l ;//(pts[i - 1].geom_distance(pts[i], R));
   }
   if (closed and pts.size() > 1)
   {
      double l = (pts.back().geom_distance(pts.front(), R));
   // std::cout << pts.back() << "->" << pts.front() << " = " << l << std::endl;
      len += (pts.back().geom_distance(pts.front(), R));
   }
   return len;
}

// #define dist(i, j) sqrt(path[i].squared_distance(path[j]))
#define dist(i, j, R) path[i].geom_distance(path[j], R)
/// - function -----------------------------------------------------------------
void two_opt(CoordsVector &path, std::vector<int> &sequence, const double R)
{
   const int N = path.size();
   int counter = 0;
   double mchange;
   do
   {
      mchange = 0.0;
      int mi = -1;
      int mj = -1;
      for (int i = 1; i < N; ++i)
      {
         for (int j = i + 1; j < N - 1; ++j)
         {
            // double change=   dist(i - 1, j) + dist(i, j + 1) - dist(i - 1, i) - dist(j, j + 1);
            double change=   dist(i - 1, j, R) + dist(i, j + 1, R) - dist(i - 1, i, R) - dist(j, j + 1, R);
         
            if (mchange > change)
            {
               mchange = change;
               mi = i;
               mj = j;
            }
            counter += 1;
         }
      }
      if (mi > 0 and mj > 0)
      {
         CoordsVector newPath;
         std::vector<int> newSequence;
         for (int i = 0; i < mi; ++i)
         {
            newPath.push_back(path[i]);
            newSequence.push_back(sequence[i]);
         }
         for (int i = mj; i >= mi; --i)
         {
            newPath.push_back(path[i]);
            newSequence.push_back(sequence[i]);
         }
         for (int i = mj + 1; i < N; ++i)
         {
            newPath.push_back(path[i]);
            newSequence.push_back(sequence[i]);
         }
         path = newPath;
         sequence = newSequence;
      }
   } while (mchange < -1e-5);
}

/* end of route_path_utils.cc */
