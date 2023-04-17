/// This is the "Cook's membrane" benchmark solved using the nonlinear elasticity solver.
/// The problem description and reference solutions can be found in the Ph.D. thesis of O.Weeger
/// "Isogeometric Finite Element Analysis of Nonlinear Structural Vibrations", 2015.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler_elasticSurface.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsIterative_multiStep.h>
#include <gsElasticity/gsWriteParaviewMultiPhysicsExtension.h>
#include <gsNurbs/gsNurbsCreator.h>

using namespace gismo;


int main(int argc, char* argv[])
{
	
	//geometry
	real_t ch = 26.0;
	real_t cb = 30.0;
	real_t bh = 30;
	real_t bb = bh;

	//material
	real_t youngsModulus = 3e-9;
	real_t poissonsRatio = 0.495;

	real_t surfaceTension = 1e-7;
	real_t surfaceYoungsModulus = 0.0;
	real_t surfacePoissonsRatio = 0.495;

	//mesh
	index_t numUniRef = 2;
	index_t numDegElev = 1;
	
	//solver
	index_t numLoadSteps = 50;
	index_t numFrames = numLoadSteps;
	index_t maxIter = 15;

	//plot
	bool meshPlot = true;
	index_t numPlotPoints = 1000;

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
	try { cmd.getValues(argc, argv); }
	catch (int rv) { return rv; }

	gsKnotVector<> KV1(0, 1, 0, 2);
	gsKnotVector<> KV2(0, 1, 0, 2);
	gsKnotVector<> KV3(0, 1, 0, 2);
	
	gsMatrix<real_t> C(8, 3);
	std::vector<gsTensorBSpline<3, real_t>> patches;
	//1
	gsVector<real_t, 3> center;
	real_t hafbx = cb / 2.;
	real_t hafby = cb / 2.;
	real_t hafh = ch / 2.;
	center << 0., 0., hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby,  hafh, hafbx, -hafby,  hafh, -hafbx, hafby,  hafh, hafbx, hafby,  hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch1(KV1, KV2, KV3, C);
	patches.push_back(patch1);
	
	//2
	hafh = bh / 2.;
	center << 0., 0., -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch2(KV1, KV2, KV3, C);
	patches.push_back(patch2);

	//3
	hafbx = bb / 2.;
	hafby = bb / 2.;
	hafh = bh / 2.;
	center << -(cb / 2. + hafbx), -(cb / 2. + hafby), -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch3(KV1, KV2, KV3, C);
	patches.push_back(patch3);

	center << (cb / 2. + hafbx), -(cb / 2. + hafby), -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch4(KV1, KV2, KV3, C);
	patches.push_back(patch4);

	center << (cb / 2. + hafbx), (cb / 2. + hafby), -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch5(KV1, KV2, KV3, C);
	patches.push_back(patch5);

	center << -(cb / 2. + hafbx), (cb / 2. + hafby), -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch6(KV1, KV2, KV3, C);
	patches.push_back(patch6);

	//4
	hafbx = cb / 2.;
	center << 0.0, -(cb / 2. + hafby), -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch7(KV1, KV2, KV3, C);
	patches.push_back(patch7);

	center << 0.0, (cb / 2. + hafby), -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch8(KV1, KV2, KV3, C);
	patches.push_back(patch8);

	//5
	hafbx = bb / 2.;
	hafby = cb / 2.;
	center << -(cb / 2. + hafbx), 0.0, -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch9(KV1, KV2, KV3, C);
	patches.push_back(patch9);

	center << (cb / 2. + hafbx), 0.0, -hafh;
	C <<
		-hafbx, -hafby, -hafh, hafbx, -hafby, -hafh, -hafbx, hafby, -hafh, hafbx, hafby, -hafh,
		-hafbx, -hafby, hafh, hafbx, -hafby, hafh, -hafbx, hafby, hafh, hafbx, hafby, hafh;
	C.col(0).array() += center(0);
	C.col(1).array() += center(1);
	C.col(2).array() += center(2);
	gsTensorBSpline<3, real_t> patch10(KV1, KV2, KV3, C);
	patches.push_back(patch10);

	gsMultiPatch<> geometry;
	for (int i = 0; i < patches.size(); i++)
		geometry.addPatch(patches[i]);

	geometry.computeTopology();

	// creating bases
	gsMultiBasis<> basisDisplacement(geometry);

	for (index_t i = 0; i < numDegElev; ++i)
		basisDisplacement.degreeElevate();
	for (index_t i = 0; i < numUniRef; ++i)
		basisDisplacement.uniformRefine();
	
	gsBoundaryConditions<> bcInfo;
	//for (int i = 0; i < geometry.nPatches(); i++)
	for (int i = 1; i < geometry.nPatches(); i++)
		for (int j = 0; j < 3; j++)
			bcInfo.addCondition(i, boundary::front, condition_type::dirichlet, nullptr, j);

	//for (int i = 0; i < geometry.nPatches(); i++)
	for (int i = 2; i < geometry.nPatches(); i++)
		bcInfo.addCondition(i, boundary::back, condition_type::robin, nullptr);
#if 1
	bcInfo.addCondition(0, boundary::back, condition_type::robin, nullptr);
	bcInfo.addCondition(0, boundary::west, condition_type::robin, nullptr);
	bcInfo.addCondition(2, boundary::west, condition_type::robin, nullptr);
	bcInfo.addCondition(5, boundary::west, condition_type::robin, nullptr);
	bcInfo.addCondition(8, boundary::west, condition_type::robin, nullptr);
	bcInfo.addCondition(0, boundary::east, condition_type::robin, nullptr);
	bcInfo.addCondition(3, boundary::east, condition_type::robin, nullptr);
	bcInfo.addCondition(4, boundary::east, condition_type::robin, nullptr);
	bcInfo.addCondition(9, boundary::east, condition_type::robin, nullptr);
	bcInfo.addCondition(0, boundary::north, condition_type::robin, nullptr);
	bcInfo.addCondition(4, boundary::north, condition_type::robin, nullptr);
	bcInfo.addCondition(5, boundary::north, condition_type::robin, nullptr);
	bcInfo.addCondition(7, boundary::north, condition_type::robin, nullptr);
	bcInfo.addCondition(0, boundary::south, condition_type::robin, nullptr);
	bcInfo.addCondition(2, boundary::south, condition_type::robin, nullptr);
	bcInfo.addCondition(3, boundary::south, condition_type::robin, nullptr);
	bcInfo.addCondition(6, boundary::south, condition_type::robin, nullptr);
#else
	bcInfo.addCondition(1, boundary::west, condition_type::robin, nullptr);
	bcInfo.addCondition(4, boundary::west, condition_type::robin, nullptr);
	bcInfo.addCondition(7, boundary::west, condition_type::robin, nullptr);
	bcInfo.addCondition(2, boundary::east, condition_type::robin, nullptr);
	bcInfo.addCondition(3, boundary::east, condition_type::robin, nullptr);
	bcInfo.addCondition(8, boundary::east, condition_type::robin, nullptr);
	bcInfo.addCondition(3, boundary::north, condition_type::robin, nullptr);
	bcInfo.addCondition(4, boundary::north, condition_type::robin, nullptr);
	bcInfo.addCondition(6, boundary::north, condition_type::robin, nullptr);
	bcInfo.addCondition(1, boundary::south, condition_type::robin, nullptr);
	bcInfo.addCondition(2, boundary::south, condition_type::robin, nullptr);
	bcInfo.addCondition(5, boundary::south, condition_type::robin, nullptr);
#endif

	//bcInfo.addCondition(2, boundary::south, condition_type::robin, nullptr);
	//bcInfo.addCondition(2, boundary::east, condition_type::robin, nullptr);
	//bcInfo.addCondition(2, boundary::west, condition_type::robin, nullptr);

	//bcInfo.addCondition(1, boundary::north, condition_type::robin, nullptr);
	//bcInfo.addCondition(1, boundary::east, condition_type::robin, nullptr);
	//bcInfo.addCondition(1, boundary::west, condition_type::robin, nullptr);


	// source function, rhs
	gsConstantFunction<> g(0., 0., 0., 3);
	
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
	gsField<> displacementField(assembler.patches(), displacement);
	gsPiecewiseFunction<> stresses;
	assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
	gsField<> stressField(assembler.patches(), stresses, true);
	
	// creating a container to plot all fields to one Paraview file
	std::map<std::string, const gsField<>*> fields;
	std::map<std::string, const gsField<>*> fields2;
	fields["Displacement"] = &displacementField;
	fields["von Mises"] = &stressField;

	// paraview collection of time steps
	std::string filenameParaview = "SurfaceSmoothing_" + util::to_string(numUniRef) + "_" + util::to_string(youngsModulus) + "_" + util::to_string(poissonsRatio) + "_" ;
	gsParaviewCollection collection(filenameParaview);
	
	std::ofstream of;
	if (numPlotPoints > 0)
	{
		gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, 0, numPlotPoints, meshPlot);
	}
	
#if 1
	gsInfo << "Solving...\n";
	index_t numStepsPerFrame = numLoadSteps / numFrames;
	index_t cs = 0;
	index_t frame = 0;
	gsStopwatch clock;
	clock.restart();
	
	for (int i = 0; i < numLoadSteps; i++)
	{
		gsInfo << "Step " << i + 1 << " of " << numLoadSteps;
		
		assembler.options().setReal("SurfaceTension", (i + 1) * surfaceTension / numLoadSteps);
		gsInfo << " with SurfaceTension: " << (i + 1) * surfaceTension / numLoadSteps << "\n";
			
		assembler.refresh();
		solver.solve(numLoadSteps);

#if 0
		gsSparseSolver<>::SimplicialLDLT LDLTsolver(assembler.matrix());
		gsVector<> stabilityVec = LDLTsolver.vectorD();
		real_t stability = stabilityVec.colwise().minCoeff()[0];
		gsInfo << "Stability = " << stability << "\n";
		if (stability < 0)
		{
			gsInfo << "Bifurcation point found!\n";
			break;
		}
#endif

		assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
		assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
		cs++;
		if (numPlotPoints > 0 && cs == numStepsPerFrame)
		{
			frame++;
			gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, frame, numPlotPoints, meshPlot);
			cs = 0;
		}
	}
	gsInfo << "Solved the system in " << clock.stop() << "s.\n";
#endif
	//delete basis;
	if (numPlotPoints > 0)
	{
		collection.save();
	}
	return 0;
}
