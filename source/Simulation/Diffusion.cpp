#include <Simulation/Diffusion.h>

void crankNicholson (Ref<VectorXd> u, const double sigma) {
    int N = u.rows();
    BandMatrix<double, N, 1, 1> A;
}
