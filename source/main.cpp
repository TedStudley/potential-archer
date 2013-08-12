#include <Geometry/shapeFunction.h>
#include <Simulation/advection.h>
#include <Simulation/diffusion.h>
#include <Geometry/norms.h>
#include <Output/latex.h>

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <Eigen/Core>

#define GRID_SIZE   64
#define SIGMA       0.9
#define VELOCITY    3.0
#define END_TIME    3.0

using namespace std;
using namespace Eigen;

std::vector <void ( *) (Ref<VectorXd>) >
shapeFnVect {gaussPulse,
             squareWave,
             triangleWave,
             sineWave
            };

std::vector <std::string>
shapeFnString {"Gaussian Pulse",
               "Square Wave",
               "Triangle Wave",
               "Sine wave"
              };

std::vector <void ( *) (Ref<VectorXd>, const double) >
methodFnVect {upwind,
              frommMethod,
              frommVanLeer
             };

std::vector <std::string>
methodFnString {"Upwind Method",
                "Fromm's Method",
                "Fromm Van Leer Method"
               };
int main() {
  double  sigma     = SIGMA,
          v         = VELOCITY,
          end_time  = END_TIME;
  std::ofstream *docStream = startTexDoc ("Tables.tex");
  for (unsigned int shapeN = 0; shapeN < shapeFnString.size(); ++shapeN) {
    (*docStream) << "\t\\section{" << shapeFnString[shapeN] << "}" << std::endl;
    for (unsigned int methodN = 0; methodN < methodFnString.size(); ++methodN) {
      (*docStream) << "\t\\subsection{" << methodFnString[methodN] << "}" << std::endl;
      makeTexTable (docStream,
                    sigma, v, end_time,
                    shapeFnVect[shapeN],
                    methodFnVect[methodN]);
    }
  }
  endTexDoc (docStream);
  return 0;
}
