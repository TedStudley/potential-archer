#pragma once

#include <cmath>
#include <iostream>
#include <Eigen/Core>

using namespace Eigen;

inline void squareWave (Ref<VectorXd> tempVect) {
  int n    = tempVect.rows();
  double h = 1.0 / n,
         x = 0.5 * h;
  for (int i = 0; i < n; ++i) {
    tempVect[i] = (0.25 < x && x < 0.75) ? 1 : 0;
    x += h;
  }
}

inline void triangleWave (Ref<VectorXd> tempVect) {
  int n    = tempVect.rows();
  double h = 1.0 / n,
         x = 0.5 * h,
         a, b;
  for (int i = 0; i < n; ++i) {
    if (x < 0.25) {
      a = 4;
      b = 0;
    } else if (x < 0.75) {
      a = -4;
      b = 2;
    } else {
      a = 4;
      b = -4;
    }
    tempVect[i] = x * a + b;
  }
}

inline void sineWave (Ref<VectorXd> tempVect) {
  int n = tempVect.rows();
  double h = 1.0 / n,
         x = 0.5 * h;
  for (int i = 0; i < n; ++i) {
    tempVect[i] = sin (x * 2 * M_PI);
    x += h;
  }
}

inline void gaussPulse (Ref<VectorXd> tempVect) {
  int n = tempVect.rows();
  double h = 1.0 / n,
         x = 0.5 * h;
  for (int i = 0; i < n; ++i) {
    tempVect[i] = exp (-256 * (x - 0.5) * (x - 0.5));
    x += h;
  }
}
