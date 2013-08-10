#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>

#include <Geometry/norms.h>
#include <Geometry/shapeFunction.h>
#include <Simulation/advection.h>

#include <vector>

std::ofstream * startTexDoc (const std::string);

void makeTexTable (std::ofstream *, 
                   const double, 
                   const double, 
                   const double,
                   void (*) (Ref<VectorXd>),
                   void(*)(Ref<VectorXd>, const double));

void endTexDoc (std::ofstream *);

