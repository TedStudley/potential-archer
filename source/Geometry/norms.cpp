#include <Geometry/norms.h>

double maxNorm (VectorXd vect) {
  return vect.maxCoeff();
}

double oneNorm (VectorXd vect) {
  return (vect.sum() / vect.rows());
}

double twoNorm (VectorXd vect) {
  return vect.norm();
}
