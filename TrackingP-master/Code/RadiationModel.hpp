#ifndef RADIATION_MODEL_HPP
#define RADIATION_MODEL_HPP

#include "Sundance.hpp"


/** --------- Model for radiation transport. Set radiation model parameters here ---- */
class RadiationModel
{
public:
  /** K0 is the light intensity at the surface */
  double K0 = 2.0;
  
  /** kappa() computes the absorption coefficient given v and u */
  Expr kappa(const Expr& pops) const ;
  
  /** kappa_u is the absorption coefficient algal biomass */
  double kappa_u = 0.0004;

  /** kappa_bg is the background turbidity due to nonphytoplankton components */
  double kappa_bg = 0.3;
  
  /** sigma_a is the absorption cross section for an algae particle (in our units) */ 
  // double sigma_a = 0.0; //0.04;

  /** Register model parameters with the command-line processor */ 
  void registerParameters();
  
  /** */
  Expr initialConditions(const Expr& popStart) const ;
};

#endif
