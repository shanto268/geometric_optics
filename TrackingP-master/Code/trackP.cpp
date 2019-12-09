#include "RadiationModel.hpp"
#include "StoichModel.hpp"
#include "TransportModel.hpp"
#include "OutputManager.hpp"
#include "DiscretizationControl.hpp"

/*
 * ============================================================
 *           Main program 
 * ============================================================
 */

int main(int argc, char** argv)
{
  try
  {
    /* ---------- Create the objects that define model parameters, etc ----------------*/
    DiscretizationControl disc;
    TransportModel trans;
    RadiationModel rad;
    LKEModel lke;
    OutputManager out;
    
    /* Register parameters to be set with command line flags */
    disc.registerParameters();
    trans.registerParameters();
    rad.registerParameters();
    lke.registerParameters();
    out.registerParameters();

    /* System initialization and command line parsing */
    Sundance::init(&argc, &argv);

    /* Check on number of populations */
    TEUCHOS_TEST_FOR_EXCEPT(lke.nPops < 1 || lke.nPops > 2);

    /* Initialize time history output */
    bool writeTimeHist=false;
    if (disc.nx<=2) writeTimeHist = true;
    ofstream of((out.dir + "/" + out.filename + "-timehist.dsv").c_str());
    if (writeTimeHist) Out::root() << "writing time history to " << out.dir << "/"
				   << out.filename << "-timehist.dsv" << endl;

    /* ------------------ We will do our linear algebra using Epetra ----------------- */
    VectorType<double> vecType = new EpetraVectorType();

    /* --------------------------- Define the problem's geometry  -------------------- */
    
    /* Get the mesh */
    Mesh mesh = disc.getMesh();

    /* Create a cell filter that will identify the maximal cells
     * in the interior of the domain */
    CellFilter interior = new MaximalCellFilter();
    /* Create a cell filter that will identify the cells at the top of the domain */
    CellFilter sides = new DimensionalCellFilter(mesh.spatialDim()-1);
    CellFilter top = sides.coordSubset(0, 0.0);
      
    /* --------------------------- Define test and unknown functions -------------------- */
    
    /* Use piecewise linear Lagrange basis for all functions */
    BasisFamily bas = new Lagrange(1);

    /* Trial and test functions for plankton */
    Expr popUnks = lke.unks();
    Expr popTests = lke.tests();
    Expr v = popUnks[0];
    Expr u = popUnks[1];
    Expr vTest = popTests[0];
    Expr uTest = popTests[1];

    /* Trial and test functions for radiation intensity */
    Expr K = new UnknownFunction(bas, "K");
    Expr KTest = new TestFunction(bas, "KTest");

    /* Trial and Test functions for algae phosphorus density */
    Expr Pa = new UnknownFunction(bas, "Pa");
    Expr PaTest = new TestFunction(bas, "PaTest");
    
    /* Trial and Test functions for free media phosphorus density */
    Expr Pf = new UnknownFunction(bas, "Pf");
    Expr PfTest = new TestFunction(bas, "PfTest");

    
    /* Aggregate functions into flat lists of expressions */
    Expr tests = List(popTests, KTest, PaTest, PfTest).flatten();
    Expr unks = List(popUnks, K, Pa, Pf).flatten();

    /* -------------------------- Define some useful symbolic expressions ------------- */
    
    /* Create differential operator and coordinate function */
    Expr grad = gradient(1);
    Expr y = new CoordExpr(0);
    Expr dy = new Derivative(0);

    /* -------------------------- Set up initial conditions --------------------------------------- */

    /* Get expressions for the ICs */
    Expr popStart = lke.initialConditions();
    Expr radStart = rad.initialConditions(popStart);
    Expr PaStart = 0.0165;
    Expr PfStart = 0.017;
    
    /* Set up the discrete space */
    Array<BasisFamily> basArray(4+lke.nPops);
    for (int i=0; i<basArray.size(); i++) basArray[i] = bas;
    DiscreteSpace discSpace(mesh, basArray, vecType);

    /* Project ICs onto the discrete space */
    L2Projector projector(discSpace,
			  List(popStart, radStart, PaStart, PfStart).flatten());
    Expr unkPrev = projector.project();

    /* Pull individual functions out of the flat list of functions */
    Expr vPrev = unkPrev[0];
    Expr uPrev = List(unkPrev[1]);
    for (int i=1; i<lke.nPops; i++) uPrev.append(unkPrev[i+1]);
    Expr KPrev = unkPrev[lke.nPops + 1];
    Expr PaPrev = unkPrev[lke.nPops + 2];
    Expr PfPrev = unkPrev[lke.nPops + 3];
    
    /* We need another discrete function for the current Newton approximation */
    Expr unkNewt = copyDiscreteFunction(unkPrev, "unkNewt");


    /* ------------------- Projector for food quality ----------------------- */
    L2Projector qProj(DiscreteSpace(mesh, new Lagrange(1), vecType),
		      lke.Q(vPrev, uPrev, PaPrev, PfPrev));
    

    /* ------------------- Set various time-related variables ----------------------- */
    
    int nSteps = disc.nSteps;
    double tFin = disc.tFin;
    double dt = tFin/((double) nSteps);
    /* Represent the time variable as a parameter expression, NOT as
     * a double variable. The reason is that we need to be able to update
     * the time value without rebuilding expressions. */
    Expr t = new Sundance::Parameter(0.0);
    Expr tPrev = new Sundance::Parameter(0.0);

    /* --------------- Define the weak form, semidiscretized in time -------------*/
    
    /* We need a quadrature rule for doing the integrations */
    QuadratureFamily quad = new GaussianQuadrature(4);
    
    /* Equation for algae dynamics */
    Expr Ra = lke.Ra(t, v, u, K, Pa, Pf);
    Expr RaPrev = lke.Ra(tPrev, vPrev, uPrev, KPrev, PaPrev, PfPrev);
    Expr algaeEqn = Integral(interior, vTest*(v - vPrev) 
			+ dt/2.0*(
				  trans.algaeTransport(vTest, v, u, Pa)
				  + trans.algaeTransport(vTest, vPrev, uPrev, PaPrev)
				  - vTest*(Ra + RaPrev)), quad);

    /* Equation for Daphnia dynamics */
    Expr Rd = lke.Rd(t, v, u, K, Pa, Pf) ;
    Expr RdPrev = lke.Rd(tPrev, vPrev, uPrev, KPrev, PaPrev, PfPrev);
    Expr Td = trans.daphniaTransport(uTest, v, u, Pa);
    Out::root() << "uPrev=" << uPrev << endl;
    Expr TdPrev = trans.daphniaTransport(uTest, vPrev, uPrev, PaPrev);
    Expr daphniaEqn = 0.0;
    
    for (int i=0; i<lke.nPops; i++)
      {
        daphniaEqn= daphniaEqn + Integral(interior, uTest[i]*(u[i]-uPrev[i]) 
			  + dt/2.0*(Td[i] + TdPrev[i]
				    - uTest[i]*(Rd[i] + RdPrev[i])
				    ), quad);
      }
    

    /* Equation for algae phosphorus transport */
    Expr RPa = lke.RPa(t, v, u, K, Pa, Pf);
    Expr RPaPrev = lke.RPa(tPrev, vPrev, uPrev, KPrev, PaPrev, PfPrev);
    Expr PaTransEqn = Integral(interior, PaTest*(Pa - PaPrev) 
			+ dt/2.0*(
				  trans.PaTransport(PaTest, v, u, Pa, Pf)
				  + trans.PaTransport(PaTest, vPrev, uPrev, PaPrev, PfPrev)
				  - PaTest*(RPa + RPaPrev)), quad);
    
    /* Equation for free media phosphorus transport */
    Expr RPf = lke.RPf(t, v, u, K, Pa, Pf);
    Expr RPfPrev = lke.RPf(tPrev, vPrev, uPrev, KPrev, PaPrev, PfPrev);
    Expr PfTransEqn = Integral(interior, PfTest*(Pf - PfPrev) 
			+ dt/2.0*(
				  trans.PfTransport(PfTest, v, u, Pa, Pf)
				  + trans.PfTransport(PfTest, vPrev, uPrev, PaPrev, PfPrev)
				  - PfTest*(RPf + RPfPrev)), quad);
    
    /* Equation for radiation transport */
    Expr radEqn= Integral(interior, KTest*(dy*K + rad.kappa(popUnks)*K), quad);

    /* Combine the five weak forms into one */
    Expr eqn = algaeEqn + daphniaEqn + radEqn + PaTransEqn + PfTransEqn;
    
    /* BC: set light intensity at top to K0. That is, K(top)-K0 = 0. 
     * The other variables,
     * v, u, Pa and Pf have no-flux ("natural") BCs and are therefore implicit in the 
     * weak form */
    Expr bc = EssentialBC(top, KTest*(K-rad.K0), quad);

    /* ------------------------ Set up nonlinear problem ---------------------------------------- */
    
    /* Here's the NLP that must be solved at each timestep. The variable uNewt 
     * is the initial guess, and when the solve is done the solution will be 
     * written into uNewt */
    NonlinearProblem prob(mesh, eqn, bc, tests, unks, unkNewt, vecType); 

    /* Choose a solver. The file "playa-newton-amesos.xml" specifies 
     * Newton-Armijo for the nonlinear solver, with Amesos (a sparse LU package)
     * for solving the linear subproblems. No need to change this until we look at
     * much larger problem sizes. */
    NonlinearSolver<double> solver 
      = NonlinearSolverBuilder::createSolver("playa-newton-amesos.xml");

    /* ------------------- Solve the problem! -----------------------------------------------------*/

    /* Write the initial conditions */
    Expr Q0 = qProj.project();
    /*Expr freeP = pfProj.project();*/
    if (writeTimeHist)
      writeTimeHistory(of, 0.0, unkPrev, Q0);
    else
      out.write(0, mesh, unkPrev, Q0);


    /* main loop over timesteps */
    Out::root() << "timestep #";
    for (int i=0; i<nSteps; i++)
    {
      /* Set the times t_i and t_{i+1} */
      if (i%100==0) Out::root() << "[" << i << "]" << endl;
      t.setParameterValue((i+1)*dt);
      tPrev.setParameterValue(i*dt);

      /* Solve the semidiscretized equations to get the solution 
       * value at the next timestep */
      SolverState<double> state = prob.solve(solver);

      /* Bail out if the solver failed */
      TEUCHOS_TEST_FOR_EXCEPTION(state.finalState() != SolveConverged,
        std::runtime_error,
        "Nonlinear solve failed to converge: message=" << state.finalMsg());
      
      Expr TotalP = Integral(interior, lke.IntegrandPt(vPrev, uPrev, PaPrev, PfPrev), quad);
      double TotalPVal = evaluateIntegral(mesh,TotalP);
     
       Out::root()<<TotalPVal<< endl;
      
      /* Set uPrev to the current solution value */
      updateDiscreteFunction(unkNewt, unkPrev);

      /* */
      Q0 = qProj.project();
      
      /* Write the current solution */
      if (writeTimeHist)
	writeTimeHistory(of, t.getParameterValue(), unkPrev, Q0);
      else if ((i+1)%disc.snapInterval==0)
	out.write((i+1)/disc.snapInterval, mesh, unkPrev, Q0);      
    } /* end of timestep loop */
  }
  catch(std::exception& e) /* Deal with errors in a reasonable graceful way */
  {
    Sundance::handleException(e);
    return -1;
  }
  Sundance::finalize(); /* system cleanup steps */
  return 0;
}

/*
 * ============================================================
 *          End of main program 
 * ============================================================
 */
