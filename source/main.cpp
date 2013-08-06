#include <Main.h>

#include <Geometry/ShapeFunction.h>
#include <Simulation/Advection.h>

int main() {
	int n = GRID_SIZE;
	double h 		= 1.0 / n,
		   sigma 	= SIGMA,
		   v		= VELOCITY,
           k        = sigma * h / v;
    int end_step    = END_TIME / k;
	Eigen::VectorXd heatVect(n);

	squareWave(heatVect);

	Eigen::VectorXd exactSoln = heatVect;

    frommVanLeer(heatVect, sigma, end_step);

	std::cout << heatVect << std::endl;

	return 0;
}
