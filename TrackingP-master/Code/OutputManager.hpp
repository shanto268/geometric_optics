#ifndef OUTPUT_MANAGER_HPP
#define OUTPUT_MANAGER_HPP

#include "Sundance.hpp"
#include <fstream>
using std::ofstream;
using std::flush;


/* ----------- Object to deal with output ----------- */
class OutputManager
{
public:
  /** */
  string filename = "test"; 

  /** */
  string dir = "SpatialResults"; 

  /* Write the result of the current step */
  void write(int step, const Mesh& mesh, const Expr& u, const Expr& Q) const ;

  /** */
  void registerParameters() ;
};

void writeTimeHistory(ofstream& of, const double& t, const Expr& u,
		      const Expr& Q);


#endif
