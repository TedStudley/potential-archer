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
  double  sigma     = SIGMA;
//          v         = VELOCITY,
//          end_time  = END_TIME;
// TeX document
/*   
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
*/

// Standard solution one-timestep output
  int N = 4096;
  VectorXd tempVect (N);
  squareWave (tempVect);
  std::cout << tempVect.transpose() << std::endl;
  squareWave (tempVect);
  std::cout << tempVect.transpose() << std::endl;
  frommVanLeer (tempVect, sigma, 27 * N / sigma);
  std::cout << tempVect.transpose() << std::endl;


// Diffusion solver
/*
  int    N = 1024;
  double h = 1.0 / N,
         delta_t = 0.01 * h;;
  bool first, second, third, fourth, fifth;
  first = second = third = fourth = fifth = 1;

  VectorXd tempVect (N);
  squareWave (tempVect);

  for (int i = 0; fifth; ++i) {
    if (first) {
      first = false;
      std::cout << tempVect.transpose() << std::endl;
    } else if (second && tempVect.maxCoeff() <= 0.8) {
      second = false;
      std::cout << tempVect.transpose() << std::endl;
    } else if (third && tempVect.maxCoeff() <= 0.6) {
      third = false;
      std::cout << tempVect.transpose() << std::endl;
    } else if (fourth && tempVect.maxCoeff() <= 0.5) {
      fourth = false;
      std::cout << tempVect.transpose() << std::endl;
    } else if (fifth && tempVect.maxCoeff() <= 0.25) {
      fifth = false;
      std::cout << tempVect.transpose() << std::endl;
    }
    crankNicholson (tempVect, sigma, delta_t);
  }
*/

  return 0;
}
