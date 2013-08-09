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
          delta_t   = 0.25 * h,
          t         = 0.0;
  int     end_step  = END_TIME / delta_t;
  VectorXd vect1 (n - 1), vect2 (n - 1), vect3 (n - 1);
  VectorXd vect (n - 1);
  sineWave (vect1, 1); sineWave (vect2, 7); sineWave (vect3, 17);
  vect = vect1 + vect2 + vect3;
  cout << vect.transpose() << endl;
  for (int timestep = 0; timestep < end_step; ++timestep) {
    upwind (vect, v, h, delta_t);
    cout << vect.transpose() << endl;
    cerr << "Timestep " << timestep << ": time = " << t << endl;
    t += delta_t;
  }
  return 0;
}
