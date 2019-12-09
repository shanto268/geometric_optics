#include "StoichModel.hpp"

void LKEModel::registerParameters()
{
  Sundance::setOption("smooth", smoothingP, "Smoothing parameter");
  Sundance::setOption("np", nPops,
		      "Number of daphnia subpopulations (must be 1 or 2)");
}

Expr LKEModel::initialConditions() const
{
  Expr vStart = 0.05;
  Array<Expr> u0(nPops);
  for (int i=0; i<u0.size(); i++)
    {
      u0[i] = 0.05/nPops;
    }
  Expr uStart = new ListExpr(u0);
  return List(vStart, uStart);
}

Expr LKEModel::unks() const
{
  BasisFamily bas = new Lagrange(1);
  Expr v = new UnknownFunction(bas, "v");
  Array<Expr> u(nPops);
  for (int i=0; i<nPops; i++)
      {
	u[i] = new UnknownFunction(bas, "u[" + std::to_string(i)+ "]");
      }
  return List(v, Expr(new ListExpr(u)));
}

Expr LKEModel::tests() const
{
  BasisFamily bas = new Lagrange(1);
  Expr v = new TestFunction(bas, "vTest");
  Array<Expr> u(nPops);
  for (int i=0; i<nPops; i++)
      {
	u[i] = new TestFunction(bas, "uTest[" + std::to_string(i)+ "]");
      }
  return List(v, Expr(new ListExpr(u)));
}


Expr LKEModel::smoothMin(const Expr& x, const Expr& y) const 
{
  double p = smoothingP;
  Expr denom = pow(pow(x, p) + pow(y, p), 1.0/p);
  return x*y/denom;
}

Expr LKEModel::maxFunc(const Expr& x) const 
{
  double zer = 0.001;
  double p = smoothingP;
  Expr argum = exp(p*zer) + exp(p*x);
  return log(argum)/p;
}

Expr LKEModel::Q(const Expr& v, const Expr& u, const Expr& Pa, const Expr& Pf) const
{
  Expr rtn = Pa;
  return rtn/v;
}

Expr LKEModel::Ra(const Expr& t, const Expr& v, const Expr& u,
		  const Expr& K, const Expr& Pa, const Expr& Pf) const 
{
  Expr K1 = maxFunc(smoothMin(K, Pa/q));
  Expr rtn = b*v*(1.0 - v/K1) ;
  Expr f = dFeeding(v);
  
  for (int i=0; i<nPops; i++) rtn = rtn - c[i]*f*u[i];

  return rtn;
}

Expr LKEModel::Rd(const Expr& t, const Expr& v, const Expr& u,
		  const Expr& K, const Expr& Pa, const Expr& Pf) const 
{
  Expr f = dFeeding(v);

  Expr PaOverX = Pa/v;
  
  Array<Expr> rtn(nPops);
  for (int i=0; i<nPops; i++)
    {
      rtn[i] = eHat*maxFunc(smoothMin(1.0, PaOverX/theta[i]))*c[i]*f*u[i]
	- delta*u[i];
    }
  return new ListExpr(rtn);
}


Expr LKEModel::RPa(const Expr& t, const Expr& v, const Expr& u,
		  const Expr& K, const Expr& Pa, const Expr& Pf) const
{
  Expr g = pFeeding(Pf);
  Expr f = dFeeding(v);

  Expr rtn = cHat*g*v-dHat*Pa;
  for (int i=0; i<nPops; i++)
    {
      rtn = rtn-(Pa/v)*c[i]*f*u[i];
    }
  return rtn;  
}

Expr LKEModel:: RPf(const Expr& t, const Expr& v, const Expr& u,
	   const Expr& K, const Expr& Pa, const Expr& Pf) const
{
  Expr g = pFeeding(Pf);
  Expr f = dFeeding(v);
  Expr PaOverX = Pa/v;
  Expr rtn = -cHat*g*v+dHat*Pa;
  for(int i=0; i<nPops; i++)
    {
      rtn = rtn + theta[i]*delta*u[i]
	+(PaOverX-eHat*theta[i]*maxFunc(smoothMin(1.0, PaOverX/theta[i])))*c[i]*f*u[i]; 
    }
  return rtn;
}

Expr LKEModel:: IntegrandPt(const Expr& v, const Expr& u,
	      const Expr& Pa, const Expr& Pf) const
{
 Expr rtn = Pa + Pf;
 for (int i=0; i<nPops; i++)
   {
     rtn = rtn + theta[i]*u[i];
   }
 return rtn;  
}
