#include "RadiationModel.hpp"

void RadiationModel::registerParameters()
{
    Sundance::setOption("K0", K0, "Surface light intensity");
    Sundance::setOption("kappa_u", kappa_u, "Algal light absorbtivity");
    Sundance::setOption("kappa_bg", kappa_bg, "Background turbidity");
    // Sundance::setOption("sigma", sigma_a, "Algae cross section");
}

Expr RadiationModel::kappa(const Expr& pops) const 
{
  return  kappa_u*pops[0]*1000 + kappa_bg;
}

Expr RadiationModel::initialConditions(const Expr& popStart) const
{
  Expr y = new CoordExpr(0);
  return K0*exp(-kappa(popStart)*y); // assumes initial populations are constant
}

