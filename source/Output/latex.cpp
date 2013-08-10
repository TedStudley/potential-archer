#include <Output/latex.h>

std::ofstream * startTexDoc (const std::string fileName) {
  static std::ofstream docStream (fileName.c_str(), std::ofstream::out | std::ofstream::trunc);
  docStream << "\\documentclass[12pt]{article}" << std::endl;
  docStream << "\\begin{document}" << std::endl;
  return & docStream;
}

void makeTexTable (std::ofstream * docStream, 
                   const double sigma, 
                   const double v, 
                   const double end_time,
                   void shapeFn (Ref<VectorXd>),
                   void methodFn (Ref<VectorXd>, const double)) {
  VectorXd tempVect;
  VectorXd exactVect;
  double oldMax, oldOne, oldTwo, newMax, newOne, newTwo;
  int n, end_step;
  double h, delta_t, t;

  (*docStream) << "\t\\begin{tabular}{l|c|c|c|c|c|c}" << std::endl;
  (*docStream) << "\t\tN & Max Norm & Convergence & One Norm & Convergence & Two Norm & Convergence \\\\" << std::endl;
  (*docStream) << "\t\t& & Rate & & Rate & & Rate \\\\" << std::endl;

  oldMax = oldOne = oldTwo = 0;

  for (n = 64; n <= 16384; n *= 2) {
    h = 1.0 / n; delta_t = sigma / v * h; t = 0; end_step = end_time / delta_t;

    (*docStream) << "\t\t\\hline" << std::endl;

    tempVect.resize (n);

    shapeFn (tempVect);
    exactVect = tempVect;

    for (int timestep = 0; timestep < end_step; ++timestep) {
      methodFn (tempVect, sigma);
      t += delta_t;
    }

    newMax = maxNorm (tempVect - exactVect);
    newOne = oneNorm (tempVect - exactVect);
    newTwo = twoNorm (tempVect - exactVect);

    (*docStream) << "\t\t" << n << std::fixed <<
            " & " << std::setprecision(6) << newMax <<
            " & " << std::setprecision(2) << ((oldMax != 0) ? log (oldMax / newMax) / log (2.0) : 0) <<
            " & " << std::setprecision(6) << newOne <<
            " & " << std::setprecision(2) << ((oldOne != 0) ? log (oldOne / newOne) / log (2.0) : 0) <<
            " & " << std::setprecision(6) << newTwo <<
            " & " << std::setprecision(2) << ((oldTwo != 0) ? log (oldTwo / newTwo) / log (2.0) : 0) << 
            " \\\\" << std::endl;
    (*docStream).unsetf (std::ios_base::floatfield);
    oldMax = newMax; oldOne = newOne; oldTwo = newTwo;
  }

  (*docStream) << "\t\\end{tabular}" << std::endl;
}

void endTexDoc (std::ofstream * docStream) {
  (*docStream) << "\\end{document}" << std::endl;
  docStream->close();
}
