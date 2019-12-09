#include "TransportModel.hpp"

void TransportModel::registerParameters()
{
    Sundance::setOption("Va", vbA, "Algae buoyant velocity");
    Sundance::setOption("Vd", vbD, "Daphnia bouyant velocity");
    Sundance::setOption("Da", Da, "Algae diffusivity");
    Sundance::setOption("Dd", Dd, "Daphnia diffusivity");
    Sundance::setOption("Dp", Dd, "Phosphorus diffusivity");
    Sundance::setOption("Dz", Dz, "Effective turbulent diffusivity");
}

Expr TransportModel::algaeTransport(const Expr& test,
				    const Expr& v,
				    const Expr& u,
				    const Expr& Pa) const
{
  Expr dx = new Derivative(0);
  return (Da + Dz)*(dx*test)*(dx*v) - vbA*(dx*test)*v;
}
 
Expr TransportModel::daphniaTransport(const Expr& test,
				      const Expr& v,
				      const Expr& u,
				      const Expr& Pa) const
{
  Expr dx = new Derivative(0);
  Array<Expr> rtn(u.size());
  for (int i=0; i<u.size(); i++)
    {
      rtn[i] =  (Dd + Dz)*(dx*test[i])*(dx*u[i]) - vbD*(dx*test[i])*u[i];
    }
  return new ListExpr(rtn);
}

Expr TransportModel::PaTransport(const Expr& test,
				const Expr& v,
				const Expr& u,
				const Expr& Pa, const Expr& Pf) const
{
  Expr dx = new Derivative(0);
  return (Dz + Da)*(dx*test)*(dx*Pa);
}

Expr TransportModel::PfTransport(const Expr& test,
				const Expr& v,
				const Expr& u,
				const Expr& Pa, const Expr& Pf) const
{
  Expr dx = new Derivative(0);
  return (Dz + Dp)*(dx*test)*(dx*Pf);
}
