/// This is the "Cook's membrane" benchmark solved using the nonlinear elasticity solver.
/// The problem description and reference solutions can be found in the Ph.D. thesis of O.Weeger
/// "Isogeometric Finite Element Analysis of Nonlinear Structural Vibrations", 2015.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler_elasticSurface.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsIterative_multiStep.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is Cook's membrane benchmark with nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/cooks.xml";
    real_t youngsModulus = 2.8;
    real_t poissonsRatio = 0.4;
    real_t surfaceTension = 2.0;
    index_t numUniRef = 3;
    index_t numDegElev = 1;
    index_t numPlotPoints = 1000;
	index_t numLoadSteps = 100;
    index_t maxIter = 100;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
    cmd.addReal("p","poisson","Poisson's ratio used in the material law",poissonsRatio);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addInt("s","point","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t", "surfTen", "Constant surface tension", surfaceTension);
    cmd.addInt("n", "numSteps", "Number of load steps", numLoadSteps);
    cmd.addInt("m", "maxIter", "Max iteration number", maxIter);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basisDisplacement.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basisDisplacement.uniformRefine();

    //gsGeometry<>* pGeom = &basisDisplacement;


    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // neumann BC
    //gsConstantFunction<> f(625e4, 0., 2);
    // elasticSurface BC
    //gsConstantFunction<> surfaceTension(4., 1);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    //for (index_t d = 0; d < 2; ++d)
        //bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,d);
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    //bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);
    bcInfo.addCondition(0, boundary::east, condition_type::robin, nullptr);
    bcInfo.addCondition(0, boundary::north, condition_type::robin, nullptr);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler
    gsElasticityAssembler_elasticSurface<real_t> assembler(geometry,basisDisplacement,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio", poissonsRatio);
    //assembler.options().setReal("SurfaceTension", surfaceTension);
    assembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    //gsIterative<real_t> solver(assembler);
    gsIterative_multiStep<real_t> solver(assembler);
    solver.options().setInt("Verbosity",solver_verbosity::all);
    solver.options().setInt("Solver",linear_solver::LDLT);
    solver.options().setInt("MaxIters", maxIter);

    gsMultiPatch<> displacement;
    assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);

    // constructing an IGA field (geometry + solution)
    gsField<> displacementField(assembler.patches(), displacement);
    //gsField<> meshAndCnet(displacement, displacement);
    gsField<> stressField(assembler.patches(), stresses, true);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string, const gsField<>*> fields;
    fields["Displacement"] = &displacementField;
    fields["von Mises"] = &stressField;
    //std::map<std::string, const gsField<>*> meshAndCnetfields;
	//meshAndCnetfields["MeshAndCnet"] = &meshAndCnet;
    // paraview collection of time steps
    gsParaviewCollection collection("cooks");
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fields, "cooks", collection, 0, numPlotPoints);
        gsWriteParaviewMultiPhysics(fields, "cooks_mesh", numPlotPoints, 1, 0);
    }
	
    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    for (int i = 0; i < numLoadSteps; i++)
    {
        assembler.options().setReal("SurfaceTension", (i + 1) * surfaceTension / numLoadSteps);
        solver.solve(numLoadSteps);
        assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
        assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
        if (numPlotPoints > 0)
            gsWriteParaviewMultiPhysicsTimeStep(fields, "cooks", collection, i + 1, numPlotPoints);
    }
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";
    if (numPlotPoints > 0)
    {
        //gsWriteParaviewMultiPhysics(fields, "cooks_final", numPlotPoints);
        collection.save();
        gsInfo << "Open \"cooks.pvd\" in Paraview for visualization.\n";
    }

    

	gsMesh<> mesh;
    //pGeom->controlNet(mesh);

    // validation
    gsMatrix<> A(2,1);
    A << 1.,1.;
    A = displacement.patch(0).eval(A);
    gsInfo << "X-displacement of the top-right corner: " << A.at(0) << std::endl;
    gsInfo << "Y-displacement of the top-right corner: " << A.at(1) << std::endl;

    return 0;
}
