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

#define GRID_SIZE   2049
#define SIGMA       0.5
#define VELOCITY    3.0
#define END_TIME    3.0

using namespace std;
using namespace Eigen;

std::vector <void (*) (Ref<VectorXd>)> 
  shapeFnVect {gaussPulse, 
               squareWave, 
               triangleWave, 
               sineWave};

std::vector <std::string> 
  shapeFnString {"Gaussian Pulse", 
                 "Square Wave", 
                 "Triangle Wave", 
                 "Sine wave"};

std::vector <void (*) (Ref<VectorXd>, const double)> 
  methodFnVect{upwind,
               frommMethod,
               frommVanLeer};

std::vector <std::string> 
  methodFnString{"Upwind Method",
                 "Fromm's Method",
                 "Fromm Van Leer Method"};
int main() {
  int     n         = GRID_SIZE;
  double  h         = 1.0 / n,
          sigma     = SIGMA,
          v         = VELOCITY,
          delta_t   = SIGMA / v * h,
          t         = 0;
  int     end_step  = END_TIME / delta_t;

  std::ofstream * docStream = startTexDoc ("Tables.tex");
  for (unsigned int shapeN = 0; shapeN < shapeFnString.size(); ++shapeN) {
    (*docStream) << "\t\\section{" << shapeFnString[shapeN] << "}" << std::endl;
    for (unsigned int methodN = 0; methodN < methodFnString.size(); ++methodN) {
      (*docStream) << "\t\\subsection{" << methodFnString[methodN] << "}" << std::endl;
      makeTexTable (docStream, 
                    sigma, v, END_TIME,
                    shapeFnVect[shapeN], 
                    methodFnVect[methodN]);
    }
  }
  endTexDoc(docStream);

  return 0;
}
