#include <Main.h>
#include <Geometry/ShapeFunction.h>
#include <Simulation/Advection.h>
#include <Simulation/Diffusion.h>

#include <iomanip>

int main() {
  int     n         = GRID_SIZE;
  double  h         = 1.0 / n,
          sigma     = SIGMA,
          v         = VELOCITY,
          delta_t   = 0.25 * h,
          t         = 0.0;
  int     end_step  = END_TIME / delta_t;
  Eigen::VectorXd heatVect (n - 1);
  squareWave (heatVect);
  for (int timestep = 0; timestep < 256 ; ++timestep) {
    crankNicholson (heatVect, sigma);
    //frommVanLeer (heatVect, sigma);
    std::cout << heatVect.transpose() << std::endl;
    t += delta_t;
  }
  return 0;
}
