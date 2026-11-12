// Cardinal temperature model (CTM), scalar version
// Shared by CTM_cpp() (exported) and rates_cpp() (called directly for speed)

#ifndef ABM_CTM_H
#define ABM_CTM_H

#include <cmath>

inline double ctm1(double tt, double top, double tmin, double tmax, double yopt) {

  double y;

  if (top - tmin < (tmax - tmin) / 2) {
    y = yopt * ((tt - tmin) * std::pow(tt - tmax, 2)) /
      ((top - tmax) * ((top - tmax) * (tt - top) -
      (top - tmin) * (top + tmax - 2 * tt)));
  } else {
    y = yopt * ((tt - tmax) * std::pow(tt - tmin, 2)) /
      ((top - tmin) * ((top - tmin) * (tt - top) -
      (top - tmax) * (top + tmin - 2 * tt)));
  }

  if (y < 0) y = 0;
  if (tt <= tmin || tt >= tmax) y = 0;

  return y;
}

#endif
