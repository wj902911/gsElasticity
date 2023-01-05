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
#include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/lshape2d_3patches_thb.xml";
    index_t numRefinementLoops = 2;
    real_t youngsModulus = 2.8;
    real_t poissonsRatio = 0.4;
    real_t surfaceTension = 2;
    index_t numUniRef = 1;
    index_t numDegElev = 1;
    index_t numPlotPoints = 10000;
    index_t numLoadSteps = 100;
    index_t maxIter = 100;
	index_t numFrames = 100;
	
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
    geometry.computeTopology();
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
    //bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    bcInfo.addCondition(0, boundary::east, condition_type::robin, nullptr);
    bcInfo.addCondition(0, boundary::north, condition_type::robin, nullptr);
    bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, nullptr, 1);
    bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, nullptr, 0);
    //bcInfo.addCondition(1, boundary::north, condition_type::robin, nullptr);
    //bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, nullptr, 1);
    bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, nullptr, 0);
    bcInfo.addCondition(2, boundary::north, condition_type::robin, nullptr);
    bcInfo.addCondition(2, boundary::east, condition_type::robin, nullptr);
    // neumann BC
    //gsConstantFunction<> f(1., 0., 0., 3);
    //bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);

    // source function, rhs
    gsConstantFunction<> g(0., 0., 2);

    // --------------- set up adaptive refinement loop ---------------
    // Specify cell-marking strategy...
    MarkingStrategy adaptRefCrit = PUCA;
    // ... and parameter.
    const real_t adaptRefParam = 0.9;

    // --------------- adaptive refinement loop ---------------
    gsInfo << "Adaptive meshing...\n";
    for (int refLoop = 0; refLoop <= numRefinementLoops; refLoop++)
    {
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

        assembler.options().setReal("SurfaceTension", surfaceTension / numLoadSteps);
        solver.solve(numLoadSteps);

        gsMultiPatch<> displacement;
        assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
        gsPiecewiseFunction<> stresses;
        assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);

        //test for adp meshing
        gsExprEvaluator<> ev;
        ev.setIntegrationElements(assembler.multiBasis());
        gsExprEvaluator<>::variable is = ev.getVariable(stresses);
        ev.integralElWise(is);
        const std::vector<real_t>& eltStress = ev.elementwise();
        for (int i = 0; i < eltStress.size(); i++)
            gsInfo << eltStress[i] << " ";
        gsInfo << std::endl;

        // --------------- adaptive refinement ---------------

        //! [adaptRefinementPart]
        // Mark elements for refinement, based on the computed local errors and
        // the refinement-criterion and -parameter.
        std::vector<bool> elMarked;
        gsMarkElementsForRef(eltStress, adaptRefCrit, adaptRefParam, elMarked);
        for (size_t k = 0; k != elMarked.size(); k++)  
            gsInfo << " " << elMarked[k];
        gsInfo << "\n";

        // Refine the marked elements with a 1-ring of cells around marked elements
        gsRefineMarkedElements(basisDisplacement, elMarked, 1);
        //! [adaptRefinementPart]

        //! [repairInterfaces]
        // Call repair interfaces to make sure that the new meshes
        // match along patch interfaces.
        basisDisplacement.repairInterfaces(geometry.interfaces());
        //! [repairInterfaces]
    }


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
    {
        gsWriteParaviewMultiPhysicsTimeStep(fields, "cube", collection, 0, numPlotPoints);
        gsWriteParaviewMultiPhysics(fields, "cube_mesh", numPlotPoints, 1, 0);
    }

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
