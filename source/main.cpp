#include <Main.h>

#include <Geometry/ShapeFunction.h>
#include <Simulation/Advection.h>
#include <Simulation/Diffusion.h>

#include <iomanip>

int main() {
	int     n = GRID_SIZE;
	double  h = 1.0 / n,
		      sigma = SIGMA,
		      v	= VELOCITY,
          k = sigma * h / v,
          t = 0.0;
  int end_step    = END_TIME / k;
	Eigen::VectorXd heatVect(n);

	squareWave(heatVect);

  std::cout << heatVect.transpose() << std::endl;
  for (int i = 0; i < 5; ++i) {
  
      crankNicholson(heatVect, sigma);
//      frommVanLeer(heatVect, sigma);
      std::cout << heatVect.transpose() << std::endl;
      t += k;
  }


	return 0;
}
