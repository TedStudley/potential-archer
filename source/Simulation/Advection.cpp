#include <Simulation/Advection.h>

// Multiple in-place upwind steps
void upwind (Ref<VectorXd> tempVect, const double sigma, const int times) {
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

// Multiple in-place Fromm's Method steps
void frommMethod (Ref<VectorXd> tempVect, const double sigma, const int times) {
	int n = tempVect.rows();
	for(int iteration = 0; iteration < times; ++iteration) {
		VectorXd oldVect = tempVect;
		for(int i = 0; i < n; ++i) {
			tempVect[i] = oldVect[i] + sigma * ((oldVect[mod (i - 1, n)] + (1.0 - sigma) * 0.25 * (oldVect[i] - oldVect[mod(i - 2, n)])) -
												(oldVect[i] + (1.0 - sigma) * 0.25 * (oldVect[mod(i + 1, n)] - oldVect[mod(i - 1, n)])));
		}
	}
}

// Multiple in-place steps of Fromm Van-Leer
void frommVanLeer (Ref<VectorXd> tempVect, const double sigma, const int times) {
	double thetaL, thetaR, phiL, phiR, VLDeltaL, VLDeltaR, fluxL, fluxR;
	int n = tempVect.rows();
	for(int iteration = 0; iteration < times; ++iteration) {
        VectorXd oldVect = tempVect;
    	for(int i = 0; i < n; ++i) {
	    	thetaL = min3 (2.0 * abs(oldVect[mod (i - 1, n)] - oldVect[mod (i - 2, n)]),
	    				   0.5 * abs(oldVect[i] - oldVect[mod (i - 2, n)]),
	    				   2.0 * abs(oldVect[i] - oldVect[mod (i - 1, n)]));
	    	thetaR = min3 (2.0 * abs(oldVect[i] - oldVect[mod (i - 1, n)]),
	    				   0.5 * abs(oldVect[mod (i + 1, n)] - oldVect[mod (i - 1, n)]),
                           2.0 * abs(oldVect[mod (i + 1, n)] - oldVect[i]));

            phiL = (oldVect[i] - oldVect[mod (i - 1, n)]) *
                   (oldVect[mod (i - 1, n)] - oldVect[mod (i - 2, n)]);
            phiR = (oldVect[mod (i + 1, n)] - oldVect[i]) *
                   (oldVect[i] - oldVect[mod (i - 1, n)]);

            VLDeltaL = (phiL > 0) ? (thetaL * sign (oldVect[i] - oldVect[mod (i - 2, n)])) : 0;
            VLDeltaR = (phiR > 0) ? (thetaR * sign (oldVect[mod (i + 1, n)] - oldVect[mod (i - 1, n)])) : 0;
        
            fluxL = oldVect[mod (i - 1, n)] + 0.5 * (1.0 - sigma) * VLDeltaL;
            fluxR = oldVect[i] + 0.5 * (1.0 - sigma) * VLDeltaR;

            tempVect[i] = oldVect[i] + sigma * (fluxL - fluxR);
        }
	}
}
