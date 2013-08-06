#pragma once

#include <cmath>
#include <Eigen/Core>
#include <Utility.h>

#define mod(x, n) 	(x + n) % n

using namespace Eigen;

// Single in-place step of upwind
inline void upwind (Ref<VectorXd> tempVect, const double sigma) {
	VectorXd oldVect = tempVect;
	double fluxL, fluxR;
	int n = tempVect.rows();
	for(int i = 0; i < n; ++i) {
		fluxL = sigma * oldVect[mod(i - 1, n)];
		fluxR = sigma * oldVect[i];
		tempVect[i] = oldVect[i] + (fluxL - fluxR);
	}
}

// Multiple in-place upwind steps
inline void upwind (Ref<VectorXd> tempVect, const double sigma, const int times) {
	double fluxL, fluxR;
	int n = tempVect.rows();
	for(int iteration = 0; iteration < times; ++iteration) {
		VectorXd oldVect = tempVect;
		for (int i = 0; i < n; ++i) {
			fluxL = sigma * oldVect[mod(i - 1, n)];
			fluxR = sigma * oldVect[i];
			tempVect[i] = oldVect[i] + (fluxL - fluxR);
		}
	}
}

// Single in-place step of Fromm's Method
inline void frommMethod (Ref<VectorXd> tempVect, const double sigma) {
	VectorXd oldVect = tempVect;
	int n = tempVect.rows();
	for(int i = 0; i < n; ++i) {
		tempVect[i] = oldVect[i] + sigma * ((oldVect[mod(i - 1, n)] + (1.0 - sigma) * 0.25 *(oldVect[i] - oldVect[mod(i - 2, n)])) -
											(oldVect[i] + (1.0 - sigma) * 0.25 * (oldVect[mod(i + 1, n)] - oldVect[mod(i - 1, n)])));
	}
}

// Multiple in-place Fromm's Method steps
inline void frommMethod (Ref<VectorXd> tempVect, const double sigma, const int times) {
	int n = tempVect.rows();
	for(int iteration = 0; iteration < times; ++iteration) {
		VectorXd oldVect = tempVect;
		for(int i = 0; i < n; ++i) {
			tempVect[i] = oldVect[i] + sigma * ((oldVect[mod (i - 1, n)] + (1.0 - sigma) * 0.25 * (oldVect[i] - oldVect[mod(i - 2, n)])) -
												(oldVect[i] + (1.0 - sigma) * 0.25 * (oldVect[mod(i + 1, n)] - oldVect[mod(i - 1, n)])));
		}
	}
}


// Single in-place step of Fromm Van-Leer
inline void frommVanLeer (Ref<VectorXd> tempVect, const double sigma) {
	double thetaL, thetaR, phiL, phiR, VLDeltaL, VLDeltaR, fluxL, fluxR;
	int n = tempVect.rows();
	VectorXd oldVect = tempVect;
	for(int i = 0; i < n; ++i) {
		thetaL = min (3,
					  2.0 * std::abs(oldVect[mod (i - 1, n)] - oldVect[mod (i - 2, n)]),
					  0.5 * std::abs(oldVect[i] - oldVect[mod (i - 2, n)]),
					  2.0 * std::abs(oldVect[i] - oldVect[mod (i - 1, n)]));
		thetaR = min (3,
					  2.0 * std::abs(oldVect[i] - oldVect[mod (i - 1, n)]),
					  0.5 * std::abs(oldVect[mod (i + 1, n)] - oldVect[mod (i - 1, n)]),
					  2.0 * std::abs(oldVect[mod (i + 1, n)] - oldVect[i]));
	}
}
