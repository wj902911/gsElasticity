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

int main(int argc, char* argv[])
{
    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/cube.xml";
    real_t youngsModulus = 2.8;
    real_t poissonsRatio = 0.4;
    real_t surfaceTension = 2.;
    index_t numUniRef = 1;
    index_t numDegElev = 1;
    index_t numPlotPoints = 1000;
    index_t numLoadSteps = 10;
    index_t maxIter = 100;
	index_t numFrames = 10;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
    cmd.addReal("p", "poisson", "Poisson's ratio used in the material law", poissonsRatio);
    cmd.addInt("r", "refine", "Number of uniform refinement application", numUniRef);
    cmd.addInt("d", "degelev", "Number of degree elevation application", numDegElev);
    cmd.addInt("s", "point", "Number of points to plot to Paraview", numPlotPoints);
    cmd.addReal("t", "surfTen", "Constant surface tension", surfaceTension);
    cmd.addInt("n", "numSteps", "Number of load steps", numLoadSteps);
    cmd.addInt("m", "maxIter", "Max iteration number", maxIter);
    cmd.addInt("f", "frames", "Number of total output frames", numFrames);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    gsMultiPatch<> geometry;
    gsReadFile<>(filename, geometry);
    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basisDisplacement.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basisDisplacement.uniformRefine();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    bcInfo.addCondition(0, boundary::front, condition_type::dirichlet, nullptr, 2);
    bcInfo.addCondition(0, boundary::east, condition_type::robin, nullptr);
    bcInfo.addCondition(0, boundary::north, condition_type::robin, nullptr);
    bcInfo.addCondition(0, boundary::back, condition_type::robin, nullptr);
    // neumann BC
    //gsConstantFunction<> f(1., 0., 0., 3);
    //bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);

    // source function, rhs
    gsConstantFunction<> g(0., 0., 0., 3);

    //=============================================//
                  // Solving & Ploting //
    //=============================================//

    // creating assembler
    gsElasticityAssembler_elasticSurface<real_t> assembler(geometry, basisDisplacement, bcInfo, g);
    assembler.options().setReal("YoungsModulus", youngsModulus);
    assembler.options().setReal("PoissonsRatio", poissonsRatio);
    assembler.options().setInt("MaterialLaw", material_law::neo_hooke_ln);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsIterative_multiStep<real_t> solver(assembler);
    solver.options().setInt("Verbosity", solver_verbosity::all);
    solver.options().setInt("Solver", linear_solver::LDLT);
    solver.options().setInt("MaxIters", maxIter);

    // constructing an IGA field (geometry + solution)
    gsMultiPatch<> displacement;
    assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
    gsField<> displacementField(assembler.patches(), displacement);
    gsField<> stressField(assembler.patches(), stresses, true);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string, const gsField<>*> fields;
    fields["Displacement"] = &displacementField;
    fields["von Mises"] = &stressField;
    // paraview collection of time steps
    gsParaviewCollection collection("cube");

    if (numPlotPoints > 0)
        gsWriteParaviewMultiPhysicsTimeStep(fields, "cube", collection, 0, numPlotPoints);
    
    gsInfo << "Solving...\n";
    index_t numStepsPerFrame = numLoadSteps / numFrames;
    index_t cs = 0;
	index_t frame = 0;
    gsStopwatch clock;
    clock.restart();
    for (int i = 0; i < numLoadSteps; i++)
    {
        assembler.options().setReal("SurfaceTension", (i + 1) * surfaceTension / numLoadSteps);
        solver.solve(numLoadSteps);
        assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
        assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
        cs++;
        if (numPlotPoints > 0 && cs == numStepsPerFrame)
        {
            frame++;
            gsWriteParaviewMultiPhysicsTimeStep(fields, "cube", collection, frame, numPlotPoints);
            cs = 0;
        }
    }
    gsInfo << "Solved the system in " << clock.stop() << "s.\n";
    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"cube.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}
