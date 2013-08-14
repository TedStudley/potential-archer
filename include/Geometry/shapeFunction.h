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

inline void fourierSquare (Ref<VectorXd> tempVect) {
  int N    = tempVect.rows(),
      m    = N / 10;
  double h = 1.0 / N,
         x = 0.5 * h;
  VectorXd bk (m);
  for (int k = 0; k < m; ++k)
    bk [k] = 2.0 * (std::cos ((k + 1) * M_PI / 4) - std::cos (3 * (k + 1) * M_PI / 4)) / ((k + 1) * M_PI);
  for (int j = 0; j < N; ++j) {
    tempVect[j] = 0;
    x = j * h + 0.5 * h;
    for (int k = 0; k < m; ++k) {
      tempVect[j] += bk[k] * std::sin ((k + 1) * M_PI * x);
    }
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
  int N = tempVect.rows() + 1;
  double h = 1.0 / N,
         x = h;
  for (int i = 1; i < N; ++i) {
    tempVect[i - 1] = sin (x * 2 * M_PI);
    x += h;
  }
}

inline void sineWave (Ref<VectorXd> tempVect, int k) {
  int N = tempVect.rows() + 1;
  double h = 1.0 / N,
         x = h;
  for (int i = 1; i < N; ++i) {
    tempVect[i - 1] = sin (x * 2 * M_PI * k);
    x += h;
  }
}

inline void gaussPulse (Ref<VectorXd> tempVect) {
  int n = tempVect.rows();
  double h = 1.0 / n,
         x = 0.5 * h;
  for (int i = 0; i < n; ++i) {
    tempVect[i] = exp (-256 * (x - 0.5) * (x - 0.5) );
    x += h;
  }
}
