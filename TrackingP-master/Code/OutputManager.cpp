#include "OutputManager.hpp"

/** */
void writeTimeHistory(ofstream& of, const double& t, const Expr& u,
		      const Expr& Q)
{
  Vector<double> vec = getDiscreteFunctionVector(u);
  Vector<double> qVec = getDiscreteFunctionVector(Q);
  of << setw(20) << t << setw(20) << vec[0] << setw(20) << vec[1]
     << setw(20) << vec[2] << setw(20) << vec[3] << setw(20) << vec[4]
     << setw(20) << qVec[1] << endl;
}

void OutputManager::registerParameters()
{
  Sundance::setOption("dir", dir, "Output directory");
  Sundance::setOption("o", filename, "Output filename");
}


void OutputManager::write(int step, const Mesh& mesh, const Expr& u,
			  const Expr& Q) const
{
  system(("mkdir -p " + dir).c_str());  

  FieldWriter writer = new DSVWriter(dir + "/" + filename + "-" 
				     + Teuchos::toString(step));
  writer.addMesh(mesh);
  writer.addField("algae", new ExprFieldWrapper(u[0]));
  for (int i=0; i<(u.size()-4); i++)
    {
      writer.addField("daphnia[" + Teuchos::toString(i)+"]",
		      new ExprFieldWrapper(u[i+1]));
    }
  writer.addField("K", new ExprFieldWrapper(u[u.size()-3]));
  writer.addField("Pa", new ExprFieldWrapper(u[u.size()-2]));
  writer.addField("Pf", new ExprFieldWrapper(u[u.size()-1]));
  writer.addField("Q", new ExprFieldWrapper(Q));
  writer.write();
}

