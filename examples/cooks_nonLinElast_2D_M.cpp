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
#include <gsElasticity/gsWriteParaviewMultiPhysicsExtension.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is Cook's membrane benchmark with nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/3patchSphereSurface.xml";
    real_t youngsModulus = 10.0;
    real_t poissonsRatio = 0.4;
    real_t surfaceTension = 1.0;
    real_t surfaceYoungsModulus = 0.0;
    real_t surfacePoissonsRatio = 0.4;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    index_t numPlotPoints = 1000;
    index_t numLoadSteps = 20;
    //index_t numPreSteps = 10;
    index_t maxIter = 100;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
    cmd.addReal("p", "poisson", "Poisson's ratio used in the material law", poissonsRatio);
    cmd.addReal("v","sufacePoisson","Poisson's ratio used in the surface material law", surfacePoissonsRatio);
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
    gsConstantFunction<> f(0., 0., 2);
	// displacement BC
    //gsConstantFunction<> f(0.1, 2);
    //gsConstantFunction<> d(-2.7e-3, 2);
    // elasticSurface BC
    //gsConstantFunction<> surfaceTension(4., 1);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    //for (index_t d = 0; d < 2; ++d)
    //{
        //bcInfo.addCornerValue(boundary::southwest, 0, 0, d);
        //bcInfo.addCornerValue(boundary::northwest, 0, 0, d);
        //bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
    //}
    //bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);
    //bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &d, 1);
    //bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 0);
    //bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    //bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);
    //bcInfo.addCondition(0, boundary::south, condition_type::robin, nullptr);
    //bcInfo.addCondition(0, boundary::east, condition_type::robin, nullptr);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler
    gsElasticityAssemblerElasticSurface<real_t> assembler(geometry,basisDisplacement,bcInfo,g);
    assembler.options().setReal("YoungsModulus",youngsModulus);
    assembler.options().setReal("PoissonsRatio", poissonsRatio);
    assembler.options().setReal("SurfaceYoungsModulus", surfaceYoungsModulus);
    assembler.options().setReal("SurfacePoissonsRatio", surfacePoissonsRatio);
    assembler.options().setReal("SurfaceTension", surfaceTension);
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
    //gsPiecewiseFunction<> stresses;
    //assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
    //gsPiecewiseFunction<> reactions;
    //assembler.constructCauchyStressesExtension(displacement, reactions, boundary::east, stress_components::all_2D_vector);
    gsPiecewiseFunction<> geoError;
    assembler.constructGeoCalc(geoError, 1.);
	
    // constructing an IGA field (geometry + solution)
    gsField<> displacementField(assembler.patches(), displacement);
    //gsField<> meshAndCnet(displacement, displacement);
    //gsField<> stressField(assembler.patches(), stresses, true);
    //gsField<> reactionField(assembler.patches(), reactions, true);
	gsField<> geoErrorField(assembler.patches(), geoError, false);
    // creating a container to plot all fields to one Paraview file
    std::map<std::string, const gsField<>*> fields;
    //std::map<std::string, const gsField<>*> fields2;
    fields["Displacement"] = &displacementField;
	fields["GeoError"] = &geoErrorField;
    //fields["von Mises"] = &stressField;
	//fields2["Reactions"] = &reactionField;
    //std::map<std::string, const gsField<>*> meshAndCnetfields;
	//meshAndCnetfields["MeshAndCnet"] = &meshAndCnet;
    // paraview collection of time steps
    std::string filenameParaview = "3patchSphereSurface_" + util::to_string(numUniRef) + "_" + util::to_string(youngsModulus) + "_" + util::to_string(surfaceYoungsModulus) + "_";
    gsParaviewCollection collection(filenameParaview);
    std::string file1 = "traction.txt";
    //std::string file2 = "displacement.txt";
    std::ofstream of;
    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, 0, numPlotPoints,true);
        gsWriteParaviewMultiPhysics(fields, "3patchSphereSurface_", numPlotPoints, false, true);
        of.open(file1);
        of << "0\n";
        of.close();

        //of.open(file2);
        //of << "0\n";
        //of.close();
		
        //gsWriteParaviewMultiPhysics(fields, filenameParaview + "mesh", numPlotPoints, 1, 0);
    }
#if 0
    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    for (int i = 0; i < numLoadSteps; i++)
    {
        //assembler.options().setReal("SurfaceTension", (i + 1) * surfaceTension / numLoadSteps);
        const gsVector<real_t> val = Eigen::Vector2d(1.61646 / numLoadSteps * (i + 1), 0.0);
		f.setValue(val, 2);
        //assembler.refresh();
        solver.solve(numLoadSteps);
        assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
        assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
        assembler.constructCauchyStressesExtension(displacement, reactions, boundary::east, stress_components::all_2D_vector);
        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, i + 1, numPlotPoints, true);
            //gsWriteHistoryOutputBoundaryResults(fields2, file1, "Reactions", boundary::east, 10);
			
			of.open(file1, std::ios_base::app);
            of << val(0) << "0\n";
            of.close();
			
            gsMatrix<> A(2, 1);
            A << 1., 1.;
            A = displacement.patch(0).eval(A);
			
            of.open(file2, std::ios_base::app);
            of << A.at(0) << "\n";
            of.close();
        }
    }
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";
#endif
    if (numPlotPoints > 0)
    {
        //gsWriteParaviewMultiPhysics(fields, "cooks_final", numPlotPoints);
        collection.save();
        gsInfo << "Open \"PlateWithHole.pvd\" in Paraview for visualization.\n";
    }

    // validation
    

    return 0;
}
