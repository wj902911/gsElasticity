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

    //std::string filename = ELAST_DATA_DIR"/3DBeam.xml";
    std::string filename = ELAST_DATA_DIR"/fitting_out.xml";
    real_t youngsModulus = 10.0;
    real_t poissonsRatio = 0.4;
    real_t surfaceTension = 0.0;
    real_t surfaceYoungsModulus = 10.0;
    real_t surfacePoissonsRatio = 0.4;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    index_t numPlotPoints = 10000;
    index_t numLoadSteps = 100;
    index_t maxIter = 100;
	index_t numFrames = 100;
	
    // minimalistic user interface for terminal
    gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
    cmd.addReal("p", "poisson", "Poisson's ratio used in the material law", poissonsRatio);
    cmd.addReal("v", "sufacePoisson", "Poisson's ratio used in the surface material law", surfacePoissonsRatio);
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
    //for (index_t i = 0; i < numDegElev; ++i)
        //basisDisplacement.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basisDisplacement.uniformRefine();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    for (index_t d = 0; d < 3; ++d)
    {
        //bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
        //bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, nullptr, d);
    }
    //bcInfo.addCondition(0, boundary::south, condition_type::robin, nullptr);
    // neumann BC
    gsConstantFunction<> f(0., -1., 0., 3);
    //bcInfo.addCondition(0, boundary::north, condition_type::neumann, &f);

    // source function, rhs
    //gsConstantFunction<> g(0., 0., 0., 3);
    gsConstantFunction<> g(0., 0., 0., 2);

    //=============================================//
                  // Solving & Ploting //
    //=============================================//

    // creating assembler
    gsElasticityAssemblerElasticSurface<real_t> assembler(geometry, basisDisplacement, bcInfo, g);
    assembler.options().setReal("YoungsModulus", youngsModulus);
    assembler.options().setReal("PoissonsRatio", poissonsRatio);
    assembler.options().setReal("SurfaceYoungsModulus", surfaceYoungsModulus);
    assembler.options().setReal("SurfacePoissonsRatio", surfacePoissonsRatio);
    assembler.options().setReal("SurfaceTension", surfaceTension);
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
    //gsPiecewiseFunction<> stresses;
    //assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
    gsField<> displacementField(assembler.patches(), displacement);
    //gsField<> stressField(assembler.patches(), stresses, true);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string, const gsField<>*> fields;
    fields["Displacement"] = &displacementField;
    //fields["von Mises"] = &stressField;
    // paraview collection of time steps
    std::string filenameParaview = "3DBeamBending_" + util::to_string(numUniRef) + "_" + util::to_string(youngsModulus) + "_" + util::to_string(surfaceYoungsModulus) + "_";
    gsParaviewCollection collection(filenameParaview);

    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fields, filenameParaview, collection, 0, numPlotPoints);
        //gsWriteParaviewMultiPhysics(fields, filenameParaview + "mesh", numPlotPoints, 1, 0);
    }
#if 0
    gsInfo << "Solving...\n";
    index_t numStepsPerFrame = numLoadSteps / numFrames;
    index_t cs = 0;
	index_t frame = 0;
    gsStopwatch clock;
    clock.restart();
    for (int i = 0; i < numLoadSteps; i++)
    {
        //assembler.options().setReal("SurfaceTension", (i + 1) * surfaceTension / numLoadSteps);
        const gsVector<real_t> val = Eigen::Vector3d(0., -1. / numLoadSteps * (i + 1), 0.);
        f.setValue(val, 3);
        solver.solve(numLoadSteps);
        assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
        assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
        cs++;
        if (numPlotPoints > 0 && cs == numStepsPerFrame)
        {
            frame++;
            gsWriteParaviewMultiPhysicsTimeStep(fields, filenameParaview, collection, frame, numPlotPoints);
            cs = 0;
        }
    }
    gsInfo << "Solved the system in " << clock.stop() << "s.\n";
    if (numPlotPoints > 0)
    {
        collection.save();
        gsInfo << "Open \"cube.pvd\" in Paraview for visualization.\n";
    }
#endif
    return 0;
}
