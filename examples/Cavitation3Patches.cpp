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

gsMatrix<real_t> generateRotationMatrix(real_t angle, gsVector<real_t> rotationAxis)
{
	gsMatrix<real_t> rotationMatrix(3, 3);
	rotationMatrix << cos(angle) + rotationAxis[0] * rotationAxis[0] * (1 - cos(angle)), rotationAxis[0] * rotationAxis[1] * (1 - cos(angle)) - rotationAxis[2] * sin(angle), rotationAxis[0] * rotationAxis[2] * (1 - cos(angle)) + rotationAxis[1] * sin(angle),
		rotationAxis[1] * rotationAxis[0] * (1 - cos(angle)) + rotationAxis[2] * sin(angle), cos(angle) + rotationAxis[1] * rotationAxis[1] * (1 - cos(angle)), rotationAxis[1] * rotationAxis[2] * (1 - cos(angle)) - rotationAxis[0] * sin(angle),
		rotationAxis[2] * rotationAxis[0] * (1 - cos(angle)) - rotationAxis[1] * sin(angle), rotationAxis[2] * rotationAxis[1] * (1 - cos(angle)) + rotationAxis[0] * sin(angle), cos(angle) + rotationAxis[2] * rotationAxis[2] * (1 - cos(angle));
	return rotationMatrix;
}

int main(int argc, char* argv[])
{
	
	//geometry
	real_t rMax = 20.0;
	real_t ir = 0.1;
	//real_t incRatial = 1.0;
	real_t x = 0.0;
	real_t y = 0.0;
	real_t z = 0.0;

	//material
	real_t youngsModulus = 1.0;
	real_t poissonsRatio = 0.495;
	real_t GA_Gr = 10.;
	
	real_t surfaceTension = GA_Gr* youngsModulus/2./(1.+ poissonsRatio)*ir;
	real_t surfaceYoungsModulus = 0.0;
	real_t surfacePoissonsRatio = 0.4;

	//mesh
	index_t numUniRef = 0;
	index_t numDegElev = 0;
	//index_t numRadialElems = 25;
	//index_t numDegRedu = 1;
	
	
	//solver
#define DISPCONTROL 1
#define FIX_INNER_SURFACE 1
#if DISPCONTROL
	real_t maxDisp = -0.07;
#else
	real_t maxLoad = 3.5;
#endif
	//index_t numPreStep = 10;
	index_t numLoadSteps = 10;
	index_t maxIter = 50;

	//plot
	bool meshPlot = true;
	index_t numFrames = 10;
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

	gsKnotVector<> KV1(0, 1, 0, 3);
	gsKnotVector<> KV2(0, 1, 0, 3);
	gsKnotVector<> KV3(0, 1, 0, 2);
	
	//gsVector<real_t> rs(numRadialElems);
	//rs.setZero();
	//rs(0) = ir * M_PI / 4. / pow(2, numUniRef) + ir;
	//for (int i = 1; i < rs.size() - 1; i++)
	    //rs(i) = rs(i - 1) * M_PI / 4. / pow(2, numUniRef) + rs(i - 1);
	std::vector<real_t> rs;
	rs.emplace_back(ir * M_PI / 4. / pow(2, numUniRef) + ir);
	rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) + rs.back());
	
#if 1
	while ((rs.back() - rs[rs.size() - 2]) * 2. + rs.back() < rMax)
	{
		if (rs.back() > (rMax - ir) / 10. + ir)
			rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) * 1.5 + rs.back());
		else if (rs.back() > (rMax - ir) / 10. * 2. + ir)
			rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) * 2.0 + rs.back());
		else
			rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) + rs.back());	
	}
#else
	while ((rs.back() - rs[rs.size() - 2]) * 2. + rs.back() < rMax)
		rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) + rs.back());
#endif
	
	rs.emplace_back(rMax);
	index_t numKnotInsert = rs.size() - 1;
	gsVector<real_t> KnotInsert;
	KnotInsert.resize(numKnotInsert);
	for (int i = 0; i < numKnotInsert; i++)
	{
		KnotInsert[i] = (rs[i] - ir) / (rMax - ir);
		KV3.insert(KnotInsert[i], 1);
	}
	//gsInfo << rs << "\n\n";
	//gsInfo << KnotInsert << "\n\n";
	gsMatrix<real_t> C(18 + 9 * numKnotInsert, 3);
	//gsMatrix<real_t> C(18, 3);
	gsMatrix<real_t> C1(9, 3);
	gsMatrix<real_t> C2(9, 3);
	gsMatrix<real_t> C3(9 * numKnotInsert, 3);
	gsMatrix<real_t> CTemp(9, 3);
	gsMatrix<real_t> W(18 + 9 * numKnotInsert, 1);
	//gsMatrix<real_t> W(18, 1);
	gsMatrix<real_t> W1(9, 1);
	gsMatrix<real_t> W2(9, 1);
	gsMatrix<real_t> W3(9 * numKnotInsert, 1);
	C3.setZero();
	W3.setZero();

	gsVector<real_t> rotAxis(3);
	rotAxis << sqrt(3.) / 3., sqrt(3.) / 3., sqrt(3.) / 3.;
	
	real_t t8 = tan(M_PI / 8.);
	real_t c8 = cos(M_PI / 8.);
	real_t c = atan(sqrt(1. / 6.) / sqrt(1. / 3.)) / 2.;
	real_t cc = cos(c);
	real_t tc = tan(c);
	real_t s12 = sqrt(1. / 2.);
	real_t s13 = sqrt(1. / 3.);
	real_t t = 1.158274940799;
	real_t X22 = 0.86385240 * t;
	real_t Y22 = 0.35620151 * t;
	real_t Z22 = 0.35620151 * t;
	real_t W22 = 0.86787082728344185507434880283629;
	
	W1 << 
		1., c8, 1.,
		c8, W22, cc,
		1., cc, 1.;
	W2 << W1;
	
	for (int i = 0; i < numKnotInsert; i++)
		W3.block(9. * i, 0, 9, 1) = W1;
	
	W << W1, W2, W3;
	//W << W1, W2;

	
	//1
	CTemp <<
		 1., 0.,  0.,  1.,  t8,  0., s12, s12,   0,
		 1., 0.,  t8, X22, Y22, Z22, s12, s12,  tc,
		s12, 0., s12, s12,  tc, s12, s13, s13, s13;
	
	C1 << CTemp * rMax;
	C2 << CTemp * ir;
	for (int i = 0; i < numKnotInsert; i++)
		C3.block(9. * i, 0, 9, 3) = CTemp * rs[i];
	//C3.block(9.*i, 0, 9, 3) = CTemp * (r* KnotInsert[i]+ ir * (1.- KnotInsert[i]));
	//C << C2, C1;
	C << C2,C3,C1;
	//gsInfo << C;

	C.col(0).array() += x;
	C.col(1).array() += y;
	C.col(2).array() += z;
	gsTensorNurbs<3, real_t> patch1(KV1, KV2, KV3, C, W);

#if 1
	//2
	gsMatrix<real_t> rotMat = generateRotationMatrix(M_PI * 2. / 3., rotAxis);
	C << C2 * rotMat.transpose(), C3* rotMat.transpose(), C1* rotMat.transpose();
	//C << C2 * rotMat.transpose(), C1* rotMat.transpose();

	gsTensorNurbs<3, real_t> patch2(KV1, KV2, KV3, C, W);
	
	//3
	rotMat = generateRotationMatrix(-M_PI * 2. / 3., rotAxis);
	C << C2 * rotMat.transpose(), C3* rotMat.transpose(), C1* rotMat.transpose();
	//C << C2 * rotMat.transpose(), C1* rotMat.transpose();


	gsTensorNurbs<3, real_t> patch3(KV1, KV2, KV3, C, W);
#endif
	gsMultiPatch<> geometry;
	geometry.addPatch(patch1);
	geometry.addPatch(patch2);
	geometry.addPatch(patch3);
	//gsInfo << patch1 << "\n\n";
	//gsInfo << patch2 << "\n\n";

	geometry.computeTopology();
	
	//gsFileData<> fd;
	//fd << geometry;
	//std::string output("geo");
	//fd.save(output);
	//gsInfo << geometry << "\n";

	// creating bases
#if 1
	gsBSplineBasis<real_t>* Bu = new gsBSplineBasis<real_t>(KV1);
	gsBSplineBasis<real_t>* Bv = new gsBSplineBasis<real_t>(KV2);
	gsBSplineBasis<real_t>* Bw = new gsBSplineBasis<real_t>(KV3);
	for (int i = 0; i < numUniRef; i++)
	{
		gsSparseMatrix<real_t, RowMajor> B[3];
		gsSparseMatrix<real_t, RowMajor> trans;
		Bu->uniformRefine_withTransfer(B[0], 1);
		Bv->uniformRefine_withTransfer(B[1], 1);
		Bw->uniformRefine_withTransfer(B[2], 0);
		tensorCombineTransferMatrices<3, real_t>(B, trans);
		W = trans * W;
	}
	gsTensorBSplineBasis<3, real_t>* tbasis = new gsTensorBSplineBasis<3, real_t>(Bu, Bv, Bw);
	gsTensorNurbsBasis<3, real_t>* basis = new gsTensorNurbsBasis<3, real_t>(tbasis, give(W));
	gsMultiBasis<> basisDisplacement;
	basisDisplacement.setTopology(geometry.topology());
	for (int i = 0; i < geometry.nPatches(); i++)
		basisDisplacement.addBasis(basis->clone());
#else
	gsMultiBasis<> basisDisplacement(geometry);
#endif


	for (index_t i = 0; i < 1; ++i)
		basisDisplacement.degreeElevate(1,2);
	for (index_t i = 0; i < numDegElev; ++i)
		basisDisplacement.degreeElevate();
	
	gsBoundaryConditions<> bcInfo;
	bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 1);
	bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 2);
	bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, nullptr, 2);
	bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, nullptr, 0);
	bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, nullptr, 0);
	bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, nullptr, 1);

	if (surfaceTension > 0)
	{
		for (int i = 0; i < geometry.nPatches(); i++)
			bcInfo.addCondition(i, boundary::front, condition_type::robin, nullptr);
	}

#if DISPCONTROL
#if FIX_INNER_SURFACE
	for (int i = 0; i < geometry.nPatches(); i++)
	{
		bcInfo.addCondition(i, boundary::front, condition_type::dirichlet, nullptr, 0);
		bcInfo.addCondition(i, boundary::front, condition_type::dirichlet, nullptr, 1);
		bcInfo.addCondition(i, boundary::front, condition_type::dirichlet, nullptr, 2);
	}
#else
	real_t dispValue = maxDisp / numLoadSteps;
	gsFunctionExpr<> dispX(util::to_string(dispValue) + "*x/" + util::to_string(ir), 3);
	gsFunctionExpr<> dispY(util::to_string(dispValue) + "*y/" + util::to_string(ir), 3);
	gsFunctionExpr<> dispZ(util::to_string(dispValue) + "*z/" + util::to_string(ir), 3);
	for (int i = 0; i < geometry.nPatches(); i++)
	{
		bcInfo.addCondition(i, boundary::front, condition_type::dirichlet, &dispX, 0);
		bcInfo.addCondition(i, boundary::front, condition_type::dirichlet, &dispY, 1);
		bcInfo.addCondition(i, boundary::front, condition_type::dirichlet, &dispZ, 2);
	}
#endif
#else
	// neumann BC
	gsFunctionExpr<> traction("10*x/0.1", "10*y/0.1", "10*z/0.1", 3);
	bcInfo.addCondition(0, boundary::front, condition_type::neumann, &traction);
	bcInfo.addCondition(1, boundary::front, condition_type::neumann, &traction);
	bcInfo.addCondition(2, boundary::front, condition_type::neumann, &traction);
#endif
	


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
	gsPiecewiseFunction<> reaction;
	assembler.constructCauchyStressesExtension(displacement, reaction, boundary::front, stress_components::normal_3D_vector);
	gsField<> reactionField(assembler.patches(), reaction, true);
	
	// creating a container to plot all fields to one Paraview file
	std::map<std::string, const gsField<>*> fields;
	std::map<std::string, const gsField<>*> fields2;
	fields["Displacement"] = &displacementField;
	fields["von Mises"] = &stressField;
	fields2["Reaction"] = &reactionField;
	//fields["CauchyStressTensor"] = &stressField;
	// paraview collection of time steps
	std::string filenameParaview = "Cavitation_" + util::to_string(numUniRef) + "_" + util::to_string(youngsModulus) + "_" + util::to_string(poissonsRatio) + "_" + util::to_string(rMax) + "_";
	gsParaviewCollection collection(filenameParaview);

#if DISPCONTROL
	std::string file1 = "disp.txt";
#else
	std::string file1 = "traction.txt";
#endif
	//std::string file2 = "InSurfEvDisp.txt";
	//std::string file2 = "InSurfEvTraction.txt";
	std::string file2 = "ouSurfEvDisp.txt";
	std::string file3 = "InSurfEvTraction.txt";
	//std::ofstream of_traction(file1);
	//std::ofstream of_inSurfEvDisp(file2);
	std::ofstream of;
	if (numPlotPoints > 0)
	{
		of.open(file1);
		of << "0\n";
		of.close();
		
		gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, 0, numPlotPoints, meshPlot);
		//gsWriteParaviewMultiPhysics(fields, filenameParaview + "mesh", numPlotPoints, 1, 0);
		
		of.open(file2);
		of << "0\n";
		of.close();
		
		of.open(file3);
		of << "0\n";
		of.close();
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
		
#if DISPCONTROL
#if FIX_INNER_SURFACE
		real_t totalDisp = 0;
		assembler.options().setReal("SurfaceTension", (i + 1) * surfaceTension / numLoadSteps);
#else
		real_t totalDisp = dispValue * (i + 1);
#endif
		of.open(file1, std::ios_base::app);
		of << totalDisp <<"\n";
		of.close();
		//dispX = gsFunctionExpr<>(util::to_string(dispValue) + "*x/" + util::to_string(ir), 3);
		//dispY = gsFunctionExpr<>(util::to_string(dispValue) + "*y/" + util::to_string(ir), 3);
		//dispZ = gsFunctionExpr<>(util::to_string(dispValue) + "*z/" + util::to_string(ir), 3);
#else
		real_t tracValue = maxLoad / numLoadSteps * (i + 1);
		of.open(file1, std::ios_base::app);
		of << tracValue << "\n";
		of.close();
		traction = gsFunctionExpr<>(
			util::to_string(tracValue) + "*x/" + util::to_string(ir),
			util::to_string(tracValue) + "*y/" + util::to_string(ir),
			util::to_string(tracValue) + "*z/" + util::to_string(ir), 3);
#endif
		
		assembler.refresh();
		solver.solve(numLoadSteps);
		assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
		assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
		assembler.constructCauchyStressesExtension(displacement, reaction, boundary::front, stress_components::normal_3D_vector);
		cs++;
		//assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
		if (numPlotPoints > 0 && cs == numStepsPerFrame)
		{
			frame++;
			gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, frame, numPlotPoints, meshPlot);
			gsWriteHistoryOutputBoundaryResults(fields, file2, "Displacement", boundary::back, 100);
			gsWriteHistoryOutputBoundaryResults(fields2, file3, "Reaction", boundary::front, 100);
			//gsWriteHistoryOutputBoundaryResults(fields, file3, "Displacement", boundary::back, numPlotPoints);
			cs = 0;
		}
	}
#endif
	gsInfo << "Solved the system in " << clock.stop() << "s.\n";
	//delete basis;
	if (numPlotPoints > 0)
	{
		collection.save();
	}
	return 0;
}
