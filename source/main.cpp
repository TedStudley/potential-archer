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
          delta_t   = 0.001 * h,
          t         = 0.0;
  int     end_step  = END_TIME / delta_t;
  firstPrint = secondPrint = thirdPrint = fourthPrint = true;
  VectorXd heatVect (n - 1);
  for (delta_t = 0.1 * h; delta_t > 0.00001 * h; delta_t *= 0.1) {
    squareWave(heatVect);
    crankNicholson (heatVect, sigma, delta_t);
    cout << heatVect.transpose() << endl;
  } 

  return 0;
}
