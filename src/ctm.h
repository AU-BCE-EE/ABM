// Cardinal temperature model (CTM), scalar version
// Used directly by rates_cpp() (for speed), and by CTM_cpp() (not exported, call as ABM:::CTM_cpp() from R,
// e.g., to plot the temperature response of a microbial group)

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
