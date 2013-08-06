#pragma once

#include <iostream>
#include <Eigen/Core>

using namespace Eigen;

inline void squareWave (Ref<VectorXd> tempVect) {
	int n = tempVect.rows();
	double h 	= 1.0 / n,
		   x	= 0.5 * h;
	for(int i = 0; i < n; ++i) {
		tempVect[i] = (0.25 < x && x < 0.75) ? 1 : 0; x += h;
	}
}
