#pragma once

#include <cmath>
#include <iostream>
#include <Eigen/Core>

using namespace Eigen;

inline void squareWave (Ref<VectorXd> tempVect) {
	int n = tempVect.rows();
	double h 	= 1.0 / n,
		   x	= 0.5 * h;
	for(int i = 0; i < n; ++i) {
		tempVect[i] = (0.485375 < x && x < 0.515625) ? 1 : 0; x += h;
	}
}

inline void sineWave (Ref<VectorXd> tempVect) {
  int n = tempVect.rows();
  double h = 1.0 / n,
         x = 0.5 * h;
  for (int i = 0; i < n; ++i) {
    tempVect[i] = sin (x * 2 * M_PI); x += h;
  }
}
