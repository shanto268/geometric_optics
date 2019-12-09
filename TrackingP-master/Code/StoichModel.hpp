#ifndef STOICH_MODEL_HPP
#define STOICH_MODEL_HPP

#include "Sundance.hpp"

/* -------------- LKE stoichiometric model ----------- */
class LKEModel
{
public:
  int nPops = 1;
  double smoothingP = 16.0;
  
  /** Food quality */
  Expr Q(const Expr& v, const Expr& u, const Expr& Pa, const Expr& Pf ) const ;

  /** */
  void registerParameters();

  /** */
  Expr Ra(const Expr& t, const Expr& v, const Expr& u,
	  const Expr& K, const Expr& Pa, const Expr& Pf) const ;

  /** */
  Expr Rd(const Expr& t, const Expr& v, const Expr& u,
	  const Expr& K, const Expr& Pa, const Expr& Pf) const ;

  /** */
  Expr RPa(const Expr& t, const Expr& v, const Expr& u,
	   const Expr& K, const Expr& Pa, const Expr& Pf) const;
  
  /** */
  Expr RPf(const Expr& t, const Expr& v, const Expr& u,
	   const Expr& K, const Expr& Pa, const Expr& Pf) const;
  
  /** */
  Expr initialConditions() const ;

  /** */
  Expr tests() const ;

  /** */
  Expr unks() const ;

  /** Integrand of Total Phosphorus*/

  Expr IntegrandPt(const Expr& v, const Expr& u,
	      const Expr& Pa, const Expr& Pf) const;

private:
  /* Smoothed min and max functions to avoid trouble during solves */
  Expr smoothMin(const Expr& x, const Expr& y) const ;
  Expr maxFunc(const Expr& x) const ;
  
  /* Feeding functions */
  Expr dFeeding(const Expr& v) const {return v/(a+v);}

  Expr pFeeding(const Expr& Pf) const {return Pf/(aHat+Pf);}
  
  
  /* LKE model parameters */

  double a = 0.25;
  double aHat = 0.008;
  double b = 1.2;
  double q = 0.0038;
  Array<double> theta = Teuchos::tuple<double>(0.03, 0.015);
  Array<double> c = Teuchos::tuple<double>(0.81, 0.5);
  double eHat = 0.8;
  double delta = 0.25;
  double cHat = 0.2;
  double dHat = 0.05;
  double T = 0.03;
};
  
#endif
