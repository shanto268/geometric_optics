#ifndef TRANSPORT_MODEL_HPP
#define TRANSPORT_MODEL_HPP

#include "Sundance.hpp"


/* --------- Model for transport. Set diffusion and buoyancy parameters here.  --------*/
class TransportModel
{
public:
  /** Dz is the turbulent diffusivity */
  double Dz = 0.1;         // m^2/day

  /** Da is the diffusivity of algae particles */
  double Da= 0.0001;     // m^2/day

  /** Dd is the diffusivity of daphnia particles */
  double Dd = 0.01;   // m^2/day

  /** Dp is the diffusivity of phosphorus particles */
  double Dp = 0.0001;   // m^2/day

  /** Buoyancy velocity of algae particles */
  double vbA = 0.0;    // meters per day

  /** Buoyancy velocity of daphnia particles */
  double vbD = 0.0;    // meters per day

  /** Register model parameters with the command-line processor */
  void registerParameters();

  /** */
  Expr algaeTransport(const Expr& test,
		      const Expr& v,
		      const Expr& u,
		      const Expr& Pa) const ;

  /** */
  Expr daphniaTransport(const Expr& test,
			const Expr& v,
			const Expr& u,
			const Expr& Pa) const ;

  /** */
  Expr PaTransport(const Expr& test,
		  const Expr& v,
		  const Expr& u,
          const Expr& Pa, const Expr& Pf) const;

  /** */
  Expr PfTransport(const Expr& test,
		  const Expr& v,
		  const Expr& u,
          const Expr& Pa, const Expr& Pf) const;
};

#endif
