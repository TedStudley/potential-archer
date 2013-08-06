#include <Main.h>

#include <Geometry/ShapeFunction.h>
#include <Simulation/Advection.h>

int main() {
	int n = GRID_SIZE;
	double h 		= 1.0 / n,
		   sigma 	= SIGMA,
		   v		= VELOCITY;
	Eigen::VectorXd heatVect(n);

	squareWave(heatVect);

	Eigen::VectorXd exactSoln = heatVect;

	frommMethod(heatVect, sigma, 1);

	std::cout << heatVect << std::endl;

	return 0;
}
