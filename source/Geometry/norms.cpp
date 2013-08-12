#include <Geometry/norms.h>

double maxNorm (VectorXd vect) {
  return vect.cwiseAbs().maxCoeff();
}

double oneNorm (VectorXd vect) {
  return (vect.cwiseAbs().sum() / vect.rows() );
}

double twoNorm (VectorXd vect) {
  return sqrt (vect.cwiseAbs2().sum() / vect.rows() );
}
