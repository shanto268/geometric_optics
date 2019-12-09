#ifndef DISCRETIZATION_CONTROL_HPP
#define DISCRETIZATION_CONTROL_HPP

#include "Sundance.hpp"


/* ----------Parameters for discretization control (space and time) --------- */
class DiscretizationControl
{
public:
  /* nx is the number of spatial mesh cells */
  int nx = 4096;
  /* nSteps is the number of timesteps to take */
  int nSteps = 1200;
  /* snapshot interval */
  int snapInterval = 5;
  /* tFin is the final simulation time. The timestep will be tFin/nSteps */
  double tFin = 240.0;
  /* Maximum depth */
  double dMax = 20.0;

  /** Register model parameters with the command-line processor */ 
  void registerParameters();

  /** Get the mesh */
  Mesh getMesh() const ;
};


#endif
