/*
 * File name: target.h
 * Date:      2016/12/09 07:55
 * Author:    Jan Faigl
 */

#ifndef __TARGET_H__
#define __TARGET_H__

#include "coords.h"

namespace gsoa {

   struct SNeuron;

   struct STarget {
      const int label;
      const Coords coords;
      const double angle;
      const double size;

      CoordsVector samples;

      int stepWinnerSelected;
      SNeuron *selectedWinner;

      STarget(const int id, const Coords &pt, const double angle, const double size) : label(id), coords(pt), stepWinnerSelected(-1), 
               selectedWinner(0), angle(angle), size(size) {}

   };

   typedef std::vector<STarget*> TargetPtrVector;

} // end namespace gsoa

#endif

/* end of target.h */
