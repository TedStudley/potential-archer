#include <Simulation/Diffusion.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

void crankNicholson (Ref<VectorXd> u, const double sigma) {
    int N = u.rows();
    
    MatrixXd A(N, N);
    A.diagonal(-1) = VectorXd::Constant(N - 1, -1);
    A.diagonal() = VectorXd::Constant(N, 2);
    A.diagonal(1) = VectorXd::Constant(N - 1, -1);
    A(0,0) = A(N-1, N-1) = 1;

    MatrixXd B = MatrixXd::Identity(N, N) + sigma * 0.5 * A;
    MatrixXd f = (MatrixXd::Identity(N, N) - sigma * 0.5 * A) * u;

    
    // Start the conjugate gradient solver

    VectorXd p, q, x, z;
    VectorXd r = f - (B * u);
    
    double rho1, rho2, beta, alpha;

    for (int i = 1;i < 30; i++) {
      z = r;
      rho1 = r.transpose() * z;
      if (i == 1) {
        p = z;
      } else {
        beta = rho1 / rho2;
        p *= beta; p += z;
      }
      q = B * p;
      alpha = rho1 / (p.transpose() * q);
      u += alpha * p;
      r -= alpha * q;
      rho2 = rho1;
      if (r.norm() < 0.00001) {
        std::cerr << "CG converged after " << i << " iterations." << std::endl;
        break;     
      }
      
   }
}

