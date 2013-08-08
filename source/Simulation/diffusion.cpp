#include <Simulation/diffusion.h>
#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Dense>

void crankNicholson (Ref<VectorXd> u, const double sigma, const double delta_t) {
  int N = u.rows() + 1;
  double h       = 1.0 / N,
         gamma = delta_t / (h * h);
  // Create the tridiagonal heat equation matrix
  MatrixXd A = MatrixXd::Zero (N - 1, N - 1);
  A.diagonal (-1) = VectorXd::Constant (N - 2, -1);
  A.diagonal() = VectorXd::Constant (N - 1, 2);
  A.diagonal (1) = VectorXd::Constant (N - 2, -1);
  A(N - 2, N - 2) = A(0, 0) = 2;
  // Set up the Crank-Nicholson equation by generating B and f
  MatrixXd B = MatrixXd::Identity (N - 1, N - 1) + (gamma * 0.5 * A);
  MatrixXd rhs = (MatrixXd::Identity (N - 1, N - 1) - gamma * 0.5 * A) * u;
  // Start the conjugate gradient solver
  VectorXd p, q, x, z;
  VectorXd r = rhs - (B * u);
  x = u;
  double rho1, rho2, beta, alpha;
  for (int i = 1;; i++) {
    z = r;
    rho1 = r.transpose() * z;
    if (i == 1) {
      p = z;
    } else {
      beta = rho1 / rho2;
      p *= beta;
      p += z;
    }
    q = B * p;
    alpha = rho1 / (p.transpose() * q);
    x += alpha * p;
    r -= alpha * q;
    rho2 = rho1;
//    std::cerr << "<==DEBUG==> r = " << r.norm() << " max = " << r.maxCoeff() << " min = " << r.minCoeff() << std::endl;
    if (r.norm() < 0.0000001) {
//      std::cerr << "CG converged after " << i << " iterations." << std::endl;
      break;
    }
  }
//  std::cerr << "<==DEBUG==> max = " << (u - x).maxCoeff() << std::endl;

  u = x;
}

