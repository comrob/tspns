/*
 * File name: coords.h
 * Date:      2013/10/13 09:23
 * Author:    Jan Faigl
 * Updated:   Jindriska Deckerova, Petr Vana
 */

#ifndef __COORDS_H__
#define __COORDS_H__

#include <vector>
#include <iostream>
#include <cmath>

/// ----------------------------------------------------------------------------
/// @brief Coords
/// ----------------------------------------------------------------------------
struct Coords
{
   double x;
   double y;
   double z;
   //double R;

   Coords(Coords &c) : x(c.x), y(c.y), z(c.z) {}
   Coords(const Coords &c) : x(c.x), y(c.y), z(c.z) {}
   Coords() {}
   Coords(double x, double y, double z) : x(x), y(y), z(z) {}

   Coords(double x, double y, double z, double R) : x(x), y(y), z(z) {}
   Coords &operator=(const Coords &c)
   {
      if (this != &c)
      {
         x = c.x;
         y = c.y;
         z = c.z;
         // R = c.R;
      }
      return *this;
   }

   inline double squared_distance(const Coords &c) const
   {
      return squared_distance(*this, c);
   }

   inline static double squared_distance(const Coords &c1, const Coords &c2)
   {
      double dx = c1.x - c2.x;
      double dy = c1.y - c2.y;
      double dz = c1.z - c2.z;
      return dx * dx + dy * dy + dz * dz;
   }

   inline Coords operator+(const Coords &c1) const
   {
      return Coords(x + c1.x, y + c1.y, z + c1.z);
   }
   inline Coords operator+(const double c) const
   {
      return Coords(x + c, y + c, z + c);
   }

   inline Coords operator-(const Coords &c1) const
   {
      return Coords(x - c1.x, y - c1.y, z - c1.z);
   }

   inline Coords operator*(const double d) const
   {
      return Coords(x * d, y * d, z * d);
   }

   inline Coords operator/(const double d) const
   {
      return Coords(x / d, y / d, z / d);
   }

   inline bool operator==(const Coords &c) const
   {
      return std::fabs(x - c.x) < 1e-03 && std::fabs(y - c.y) && std::fabs(z - c.z);
   }

   inline double dot(const Coords &c1) const
   {
      return x * c1.x + y * c1.y + z * c1.z;
   }

   inline double norm() const
   {
      return sqrt(dot(*this));
   }

   inline Coords norm2() const
   {
      return *this / sqrt(dot(*this));
   }

   inline double angle(const Coords &c1) const
   {
      double tmp = (dot(c1)) / (norm() * c1.norm());
      if (tmp >= 1)
         return 0;
      if (tmp <= -1)
         return M_PI;
      return acos(tmp);
   }

   inline double geom_distance(const Coords &c, const double R) const
   {
      return geom_distance(*this, c, R);
   }

   inline double geom_distance(const Coords &c1, const Coords &c2, const double &R) const
   {
      return R * c1.angle(c2);
   }

   inline Coords cross(Coords c1) const
   {
      return Coords(y * c1.z - z * c1.y, z * c1.x - x * c1.z, x * c1.y - y * c1.x);
   }

   inline std::vector<Coords> sample(const double &ang, const int &m, const double &R)
   {
      std::vector<Coords> out;

      Coords v1(1, 0, 0);
      if (std::fabs(dot(v1)) < 0.01)
      {
         v1 = Coords(0, 0, 1);
      }
      double sang = tan(ang);
      Coords d1 = cross(v1).norm2() * sang;
      Coords d2 = cross(d1).norm2() * sang;
      Coords x = norm2();
      for (double alpha = 0.0; alpha < 2 * M_PI; alpha += 2 * M_PI / m)
      {
         Coords p = x + d1 * cos(alpha) + d2 * sin(alpha);
         out.push_back(p.norm2() * R);
      }
      return out;
   }

   inline void generate_random()
   {
      do
      {
         x = (double)rand() / RAND_MAX;
         y = (double)rand() / RAND_MAX;
         z = (double)rand() / RAND_MAX;
      } while (norm() < 0.001);
      *this = norm2();
   }
};

// interpolate between points a,b, coef \in [0, 1]
static inline Coords interpolate(const Coords a, const Coords b, double coef)
{
   Coords x = a * (1 - coef) + b * coef;
   return x.norm2();
}

static inline Coords closest1(const Coords &p1, const Coords &p2, const Coords &target)
{
   Coords a = p1;
   if (target.angle(p1) > target.angle(p2))
   {
      a = p2;
   }
   Coords n = p1.cross(p2);
   if (n.norm() >= 0.0001)
   {
      Coords x = target.cross(n).cross(n);
      if (x.norm() > 0.0001)
      {
         if (p1.cross(n).dot(x) >= 0 && p2.cross(n).dot(x) <= 0)
         {
            a = x * -1;
         }
      }
   }
   return a;
}

static inline Coords closest2(const double &omega, const Coords &p1,
                              const Coords &p2, const Coords &target, Coords &x)
{
   x = closest1(p1, p2, target);
   if (x.angle(target) > omega)
   {
      Coords n = x.cross(target);
      double sang = tan(omega);
      Coords d1 = n.cross(target).norm2() * sang;
      if (d1.dot(x) < 0)
      {
         d1 = d1 * -1;
      }
      Coords p = target.norm2() + d1;
      return p.norm2();
   }
   return x;
}

static inline Coords closest_with_optimization(const double &omega, const Coords &p1,
                                               const Coords &p2, const Coords &target, Coords &x)
{
   x = closest1(p1, p2, target);
   if (x.angle(target) > omega)
   {
      // std::cout << "post opt " << std::endl;
      Coords n = x.cross(target);
      double sang = tan(omega); //tan(omega);
      Coords d1 = n.cross(target).norm2() * sang;
      Coords d2 = d1.cross(target).norm2() * sang;
      if (d1.dot(x) < 0)
      {
         d1 = d1 * -1;
      }

      double alpha = 0;
      double step = 0.1;

      Coords p_best = target.norm2() + d1;
      p_best = p_best.norm2();
      double l_best = p1.angle(p_best) + p2.angle(p_best);
      while (std::fabs(step) > 1e-10)
      {
         double alpha2 = alpha + step;
         Coords p = target.norm2() + d1 * cos(alpha2) + d2 * sin(alpha2);
         p = p.norm2();
         double l = p1.angle(p) + p2.angle(p);
         if (l < l_best)
         {
            l_best = l;
            step *= 2;
            p_best = p;
            alpha = alpha2;
         }
         else
         {
            step *= -0.1;
         }
      }

      return p_best;
   }
   return x;
}

typedef std::vector<Coords> CoordsVector;
typedef std::vector<CoordsVector> CoordsVectorVector;

std::ostream &operator<<(std::ostream &os, const Coords &pt);

std::istream &operator>>(std::istream &is, Coords &pt);

std::istream &operator>>(std::istream &is, CoordsVector &pts);

#endif

/* end of coords.h */
