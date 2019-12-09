#include "DiscretizationControl.hpp"

void DiscretizationControl::registerParameters()
{
    Sundance::setOption("snap", snapInterval, "Snapshot writing interval");
    Sundance::setOption("nx", nx, "Num spatial elements");
    Sundance::setOption("nt", nSteps, "Num timesteps");
    Sundance::setOption("tf", tFin, "Final time");
    Sundance::setOption("dMax", dMax, "Max depth");
}

Mesh DiscretizationControl::getMesh() const
{
  MeshType meshType = new BasicSimplicialMeshType();
  MeshSource mesher = new PartitionedLineMesher(0.0, dMax, nx, meshType);
  Mesh mesh = mesher.getMesh();
  return mesh;
}



