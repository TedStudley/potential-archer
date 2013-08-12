#pragma once

#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Utility.h>

#define mod(x, n) 	(x + n) % n
#define sign(x)     (x < 0) ? (-1) : 1
#define max2(x, y)    (x > y) ? x : y
#define max3(x, y, z) (x > y) ? ((x > z) ? x : z) : ((y > z) ? y : z)
#define min2(x, y)    (x < y) ? x : y
#define min3(x, y, z) (x < y) ? ((x < z) ? x : z) : ((y < z) ? y : z)

using namespace Eigen;

// Single in-place step of upwind
inline void upwind (Ref<VectorXd> tempVect, const double sigma) {
  VectorXd oldVect = tempVect;
  double fluxL, fluxR;
  int n = tempVect.rows();
  for (int i = 0; i < n; ++i) {
    fluxL = sigma * oldVect[mod (i - 1, n)];
    fluxR = sigma * oldVect[i];
    tempVect[i] = oldVect[i] + (fluxL - fluxR);
  }
}

// Multiple in-place upwind steps
void upwind (Ref<VectorXd>, const double, const int);

// Single in-place step of Fromm's Method
inline void frommMethod (Ref<VectorXd> tempVect, const double sigma) {
  VectorXd oldVect = tempVect;
  int n = tempVect.rows();
  for (int i = 0; i < n; ++i) {
    tempVect[i] = oldVect[i] + sigma * ( (oldVect[mod (i - 1, n)] + (1.0 - sigma) * 0.25 * (oldVect[i] - oldVect[mod (i - 2, n)]) ) -
                                         (oldVect[i] + (1.0 - sigma) * 0.25 * (oldVect[mod (i + 1, n)] - oldVect[mod (i - 1, n)]) ) );
  }
}

// Multiple in-place Fromm's Method steps
void frommMethod (Ref<VectorXd>, const double, const int);

// Single in-place step of Fromm Van-Leer
inline void frommVanLeer (Ref<VectorXd> tempVect, const double sigma) {
  double thetaL, thetaR, phiL, phiR, VLDeltaL, VLDeltaR, fluxL, fluxR;
  int n = tempVect.rows();
  VectorXd oldVect = tempVect;
  for (int i = 0; i < n; ++i) {
    thetaL = min3 (2.0 * std::abs (oldVect[mod (i - 1, n)] - oldVect[mod (i - 2, n)]),
                   0.5 * std::abs (oldVect[i] - oldVect[mod (i - 2, n)]),
                   2.0 * std::abs (oldVect[i] - oldVect[mod (i - 1, n)]) );
    thetaR = min3 (2.0 * std::abs (oldVect[i] - oldVect[mod (i - 1, n)]),
                   0.5 * std::abs (oldVect[mod (i + 1, n)] - oldVect[mod (i - 1, n)]),
                   2.0 * std::abs (oldVect[mod (i + 1, n)] - oldVect[i]) );
    phiL = (oldVect[i] - oldVect[mod (i - 1, n)]) *
           (oldVect[mod (i - 1, n)] - oldVect[mod (i - 2, n)]);
    phiR = (oldVect[mod (i + 1, n)] - oldVect[i]) *
           (oldVect[i] - oldVect[mod (i - 1, n)]);
    VLDeltaL = (phiL > 0) ? copysign (thetaL, (oldVect[i] - oldVect[mod (i - 2, n)]) ) : 0;
    VLDeltaR = (phiR > 0) ? copysign (thetaR, (oldVect[mod (i + 1, n)] - oldVect[mod (i - 1, n)]) ) : 0;
    fluxL = oldVect[mod (i - 1, n)] + 0.5 * (1.0 - sigma) * VLDeltaL;
    fluxR = oldVect[i] + 0.5 * (1.0 - sigma) * VLDeltaR;
    tempVect[i] = oldVect[i] + sigma * (fluxL - fluxR);
  }
}

// Multiple in-place steps of Fromm Van-Leer
void frommVanLeer (Ref<VectorXd>, const double, const int);
