#include <main.h>
#include <Geometry/shapeFunction.h>
#include <Simulation/advection.h>
#include <Simulation/diffusion.h>

#include <iomanip>
#include <cmath>

using namespace std;
using namespace Eigen;

int main() {
  int     n         = GRID_SIZE;
  double  h         = 1.0 / n,
          sigma     = SIGMA,
          v         = VELOCITY,
          delta_t   = 0.01 * h,
          t         = 0.0;
  int     end_step  = END_TIME / delta_t;
  VectorXd vect1 (n - 1);  // 0.00001
  sineWave (vect1);
  for (int timestep = 0; timestep < 1024; ++timestep) {
    crankNicholson (vect1, sigma, delta_t);
    cout << vect1.transpose() << endl;
  }
    return 0;
}
