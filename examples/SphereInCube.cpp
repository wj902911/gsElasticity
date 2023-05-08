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
	real_t cubeSize = 20.0;
	real_t ir = 0.1;//0.1mm
	//real_t incRatial = 1.0;
	real_t x = 0.0;
	real_t y = 0.0;
	real_t z = 0.0;

	//material
	real_t youngsModulus = 10;//10Mpa
	real_t poissonsRatio = 0.495;
//#define USE_SURFACE_TENSION
#ifdef USE_SURFACE_TENSION
	//real_t GA_Gr = 10.;
	//real_t surfaceTension = GA_Gr * youngsModulus / 2. / (1. + poissonsRatio) * ir;
	real_t surfaceTension = 0.3;
#else
	real_t surfaceTension = 0.0;
	real_t DeltaST = 0.0;
#endif // USE_SURFACE_TENSION
	real_t surfaceYoungsModulus = 0.0;
	real_t surfacePoissonsRatio = 0.0;

	//mesh
	index_t numUniRef_hoop = 1;
	index_t numUniRef_radial = 3;
	index_t numDegElev = 0;
	//index_t numRadialElems = 25;
	//index_t numDegRedu = 1;

	index_t maxBisecNum = 20;
	real_t bisecTol = 1e-10;

	//solver
#define DISPCONTROL 1
//#define FIX_INNER_SURFACE 0
#if DISPCONTROL
	real_t maxDisp = 4;
#else
	real_t maxLoad = 50.0;
#endif
#ifdef USE_SURFACE_TENSION
	index_t numPreStep = 15;
	real_t DeltaST = surfaceTension / numPreStep;
	index_t numLoadSteps = 35;
	index_t numFrames = numLoadSteps;
#else
	index_t numPreStep = 0;
	index_t numLoadSteps = 20;
	index_t numFrames = numLoadSteps;
#endif // USE_SURFACE_TENSION
#if DISPCONTROL
	real_t DeltaDisp = maxDisp / (numLoadSteps - numPreStep);
#else
	real_t DeltaTrac = maxLoad / (numLoadSteps - numPreStep);
#endif
	index_t maxIter = 30;

	//plot
	bool meshPlot = true;
	index_t numPlotPoints = 1000;

	// minimalistic user interface for terminal
	gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
	cmd.addReal("p", "poisson", "Poisson's ratio used in the material law", poissonsRatio);
	cmd.addReal("v", "sufacePoisson", "Poisson's ratio used in the surface material law", surfacePoissonsRatio);
	cmd.addInt("r", "refine_r", "Number of uniform refinement application in radial direction", numUniRef_radial);
	cmd.addInt("o", "refine_h", "Number of uniform refinement application in hoop direction", numUniRef_hoop);
	cmd.addInt("d", "degelev", "Number of degree elevation application", numDegElev);
	cmd.addInt("s", "point", "Number of points to plot to Paraview", numPlotPoints);
	cmd.addReal("t", "surfTen", "Constant surface tension", surfaceTension);
	cmd.addInt("n", "numSteps", "Number of load steps", numLoadSteps);
	cmd.addInt("m", "maxIter", "Max iteration number", maxIter);
	try { cmd.getValues(argc, argv); }
	catch (int rv) { return rv; }

	gsKnotVector<> KV1(0, 1, 0, 5);
	gsKnotVector<> KV2(0, 1, 0, 5);
	gsKnotVector<> KV3(0, 1, 0, 2);

	//gsVector<real_t> rs(numRadialElems);
	//rs.setZero();
	//rs(0) = ir * M_PI / 4. / pow(2, numUniRef) + ir;
	//for (int i = 1; i < rs.size() - 1; i++)
		//rs(i) = rs(i - 1) * M_PI / 4. / pow(2, numUniRef) + rs(i - 1);
#if 0
	real_t ratio1 = 1.5;
	real_t ratio2 = 2.0;
	real_t ratio3 = 2.0;
	std::vector<real_t> rs;
	rs.emplace_back(ir * M_PI / 4. / pow(2, numUniRef) * ratio1 + ir);
	rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) * ratio1 + rs.back());

#if 1
	while ((rs.back() - rs[rs.size() - 2]) * 2. + rs.back() < rMax)
	{
		if (rs.back() > (rMax - ir) / 10. + ir)
			rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) * ratio2 + rs.back());
		else if (rs.back() > (rMax - ir) / 10. * 2. + ir)
			rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) * ratio3 + rs.back());
		else
			rs.emplace_back(rs.back() * M_PI / 4. / pow(2, numUniRef) * ratio1 + rs.back());
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
#else
	std::vector<real_t> rs;
	index_t numKnotInsert = 0;
#endif
	//gsInfo << rs << "\n\n";
	//gsInfo << KnotInsert << "\n\n";
	gsMatrix<real_t> C(50 + 25 * numKnotInsert, 3);
	//gsMatrix<real_t> C(18, 3);
	gsMatrix<real_t> C1(25, 3);
	gsMatrix<real_t> C2(25, 3);
	//gsMatrix<real_t> C3(25 * numKnotInsert, 3);
	//gsMatrix<real_t> CTemp(25, 3);
	gsMatrix<real_t> W(50 + 25 * numKnotInsert, 1);
	//gsMatrix<real_t> W(18, 1);
	gsMatrix<real_t> W1(25, 1);
	gsMatrix<real_t> W2(25, 1);
	//gsMatrix<real_t> W3(25 * numKnotInsert, 1);
	//C3.setZero();
	//W3.setZero();

	gsVector<real_t> rotAxis(3);
	rotAxis << sqrt(3.) / 3., sqrt(3.) / 3., sqrt(3.) / 3.;

	W1 <<
		4.00000000000000000000000000000000,
		4.00000000000000000000000000000000,
		4.11438191683587284330769762163982,
		4.34314575050761941810151256504469,
		4.68629150101523972438144483021460,
		4.00000000000000000000000000000000,
		4.00000000000000000000000000000000,
		4.10867849108508043087795158498921,
		4.32456764352295053299712890293449,
		4.64619962758132043489922580192797,
		4.11438191683587284330769762163982,
		4.10867849108508043087795158498921,
		4.20526495214135032085778220789507,
		4.40174810706343944133323020651005,
		4.69731520891156240082864314899780,
		4.34314575050761941810151256504469,
		4.32456764352295053299712890293449,
		4.40174810706343944133323020651005,
		4.57237657197148017473864456405863,
		4.83730323642847892529061937239021,
		4.68629150101523972438144483021460,
		4.64619962758132043489922580192797,
		4.69731520891156240082864314899780,
		4.83730323642847892529061937239021,
		5.07179676972449122729358350625262;
	
	W2 << 
		1., 1., 1., 1., 1.,
		1., 1., 1., 1., 1.,
		1., 1., 1., 1., 1.,
		1., 1., 1., 1., 1.,
		1., 1., 1., 1., 1.;

	//for (int i = 0; i < numKnotInsert; i++)
	//	W3.block(25. * i, 0, 25, 1) = W1;

	W << W1, W2;


	//1
	C1 <<
		1.00000000000000000000000000000000, 0.00000000000000000000000000000000, 0.00000000000000000000000000000000,
		1.00000000000000000000000000000000, 0.20710678118654751722615969811159, 0.00000000000000000000000000000000,
		0.94439897941033268402577505185036, 0.40269821396808214153395510948030, 0.00000000000000000000000000000000,
		0.84198285288145657823122292029439, 0.57223070949163856724339893844444, 0.00000000000000000000000000000000,
		0.70710678118654746171500846685376, 0.70710678118654746171500846685376, 0.00000000000000000000000000000000,

		1.00000000000000000000000000000000, 0.00000000000000000000000000000000, 0.20710678118654751722615969811159,
		1.00000000000000000000000000000000, 0.20710678118654751722615969811159, 0.20710678118654751722615969811159,
		0.94432179734478771671035701729124, 0.40268222903892175734696934341628, 0.19797538630529729064555510831269,
		0.84198285288145646720892045777873, 0.57223070949163856724339893844444, 0.18115045237214510986945015247329,
		0.70710678118654746171500846685376, 0.70710678118654746171500846685376, 0.15891862259789110711771797923575,

		0.94439897941033268402577505185036, 0.00000000000000000000000000000000, 0.40269821396808214153395510948030,
		0.94432179734478771671035701729124, 0.19797538630529729064555510831269, 0.40268222903892175734696934341628,
		0.89514378285984974592537355420063, 0.38556611722930161922917591255100, 0.38556611722930161922917591255100,
		0.80450256336705749937721066089580, 0.55014710267914124219856830677600, 0.35385983263387738029237539194582,
		0.68375661096064965782659328397131, 0.68375661096064965782659328397131, 0.31168902873712167611586210114183,

		0.84198285288145657823122292029439, 0.00000000000000000000000000000000, 0.57223070949163856724339893844444,
		0.84198285288145646720892045777873, 0.18115045237214510986945015247329, 0.57223070949163856724339893844444,
		0.80450256336705749937721066089580, 0.35385983263387738029237539194582, 0.55014710267914124219856830677600,
		0.73473280857926925868639500549762, 0.50883886682844026161376405070769, 0.50883886682844026161376405070769,
		0.63966542995691610951070060764323, 0.63966542995691610951070060764323, 0.45271994765504508517750537066604,

		0.70710678118654746171500846685376, 0.00000000000000000000000000000000, 0.70710678118654746171500846685376,
		0.70710678118654746171500846685376, 0.15891862259789110711771797923575, 0.70710678118654746171500846685376,
		0.68375661096064965782659328397131, 0.31168902873712167611586210114183, 0.68375661096064965782659328397131,
		0.63966542995691610951070060764323, 0.45271994765504508517750537066604, 0.63966542995691610951070060764323,
		0.57735026918962573105886804114562, 0.57735026918962573105886804114562, 0.57735026918962573105886804114562;

	std::ofstream file("conctrolPoints.txt");
	file << std::setprecision(32);
	gsMatrix<real_t> rotMat_o = generateRotationMatrix(-M_PI * 2. / 3., rotAxis);
	file << C1 * rotMat_o.transpose();
	file.close();
	
	C1 << C1 * ir;
	C2 <<
		1., 0., 0.,
		1., .25, 0.,
		1., .5, 0.,
		1., .75, 0.,
		1., 1., 0.,

		1., 0., .25,
		1., .25, .25,
		1., .5, .25,
		1., .75, .25,
		1., 1., .25,

		1., 0., .5,
		1., .25, .5,
		1., .5, .5,
		1., .75, .5,
		1., 1., .5,

		1., 0., .75,
		1., .25, .75,
		1., .5, .75,
		1., .75, .75,
		1., 1., .75,

		1., 0., 1.,
		1., .25, 1.,
		1., .5, 1.,
		1., .75, 1.,
		1., 1., 1.;

	C2 << C2 * cubeSize;
	//for (int i = 0; i < numKnotInsert; i++)
		//C3.block(25. * i, 0, 25, 3) = CTemp * rs[i];
	//C3.block(9.*i, 0, 9, 3) = CTemp * (r* KnotInsert[i]+ ir * (1.- KnotInsert[i]));
	//C << C2, C1;
	C << C1, C2;
	//gsInfo << C;

	C.col(0).array() += x;
	C.col(1).array() += y;
	C.col(2).array() += z;
	gsTensorNurbs<3, real_t> patch1(KV1, KV2, KV3, C, W);

#if 1
	//2
	gsMatrix<real_t> rotMat = generateRotationMatrix(M_PI * 2. / 3., rotAxis);
	C << C1 * rotMat.transpose(), C2* rotMat.transpose();
	//C << C2 * rotMat.transpose(), C1* rotMat.transpose();

	gsTensorNurbs<3, real_t> patch2(KV1, KV2, KV3, C, W);

	//3
	rotMat = generateRotationMatrix(-M_PI * 2. / 3., rotAxis);
	C << C1 * rotMat.transpose(), C2* rotMat.transpose();
	//C << C2 * rotMat.transpose(), C1* rotMat.transpose();


	gsTensorNurbs<3, real_t> patch3(KV1, KV2, KV3, C, W);
#endif
	gsMultiPatch<> geometry;
	geometry.addPatch(patch1);
	geometry.addPatch(patch2);
	geometry.addPatch(patch3);

	geometry.computeTopology();

	// creating bases
#if 1
	gsBSplineBasis<real_t>* Bu = new gsBSplineBasis<real_t>(KV1);
	gsBSplineBasis<real_t>* Bv = new gsBSplineBasis<real_t>(KV2);
	gsBSplineBasis<real_t>* Bw = new gsBSplineBasis<real_t>(KV3);
#if 1
	for (int i = 0; i < numUniRef_hoop; i++)
	{
		gsSparseMatrix<real_t, RowMajor> B[3];
		gsSparseMatrix<real_t, RowMajor> trans;
		Bu->uniformRefine_withTransfer(B[0], 1);
		Bv->uniformRefine_withTransfer(B[1], 1);
		Bw->uniformRefine_withTransfer(B[2], 0);
		tensorCombineTransferMatrices<3, real_t>(B, trans);
		W = trans * W;
	}
	for (int i = 0; i < numUniRef_radial; i++)
	{
		gsSparseMatrix<real_t, RowMajor> B[3];
		gsSparseMatrix<real_t, RowMajor> trans;
		Bu->uniformRefine_withTransfer(B[0], 0);
		Bv->uniformRefine_withTransfer(B[1], 0);
		Bw->uniformRefine_withTransfer(B[2], 1);
		tensorCombineTransferMatrices<3, real_t>(B, trans);
		W = trans * W;
	}
	gsTensorBSplineBasis<3, real_t>* tbasis = new gsTensorBSplineBasis<3, real_t>(Bu, Bv, Bw);
#else
	gsTensorBSplineBasis<3, real_t>* tbasis = new gsTensorBSplineBasis<3, real_t>(Bu, Bv, Bw);
	for (index_t i = 0; i < numUniRef_hoop; i++)
	{
		tbasis->x().uniformRefine();
		tbasis->y().uniformRefine();
	}
	for (index_t i = 0; i < numUniRef_radial; i++)
		tbasis->z().uniformRefine();
#endif
	gsTensorNurbsBasis<3, real_t>* basis = new gsTensorNurbsBasis<3, real_t>(tbasis, give(W));
	gsMultiBasis<> basisDisplacement;
	basisDisplacement.setTopology(geometry.topology());
	for (int i = 0; i < geometry.nPatches(); i++)
		basisDisplacement.addBasis(basis->clone());
#else
	gsMultiBasis<> basisDisplacement(geometry);
#endif


	for (index_t i = 0; i < 2; ++i)
		basisDisplacement.degreeElevate(1, 2);
	for (index_t i = 0; i < numDegElev; ++i)
		basisDisplacement.degreeElevate();

	gsBoundaryConditions<> bcInfo;
	bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 1);
	bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 2);
	bcInfo.addCondition(0, boundary::back, condition_type::dirichlet, nullptr, 0);
	bcInfo.addCondition(1, boundary::west, condition_type::dirichlet, nullptr, 2);
	bcInfo.addCondition(1, boundary::south, condition_type::dirichlet, nullptr, 0);
	bcInfo.addCondition(1, boundary::back, condition_type::dirichlet, nullptr, 1);
	bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, nullptr, 0);
	bcInfo.addCondition(2, boundary::south, condition_type::dirichlet, nullptr, 1);
	bcInfo.addCondition(2, boundary::back, condition_type::dirichlet, nullptr, 0);
	bcInfo.addCondition(2, boundary::back, condition_type::dirichlet, nullptr, 1);

	if (surfaceTension > 0 || surfaceYoungsModulus > 0)
	{
		for (int i = 0; i < geometry.nPatches(); i++)
			bcInfo.addCondition(i, boundary::front, condition_type::robin, nullptr);
	}

#if DISPCONTROL
	gsConstantFunction<> dispZ(DeltaDisp, 3);
	for (int i = 0; i < geometry.nPatches(); i++)
	{
		bcInfo.addCondition(i, boundary::back, condition_type::dirichlet, &dispZ, 2);
	}
#else
	// neumann BC
	gsFunctionExpr<> traction("0.0", "0.0", "0.0", 3);
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
	//std::string filenameParaview = "SphereInCube_" + util::to_string(numUniRef_hoop) + "_" + util::to_string(youngsModulus) + "_" + util::to_string(poissonsRatio) + "_" + util::to_string(cubeSize) + "_";
	std::string filenameParaview = "SphereInCube";
	gsParaviewCollection collection(filenameParaview);

#if DISPCONTROL
	std::string file1 = "disp.txt";
#else
	std::string file1 = "traction.txt";
	std::string file2 = "InSurfEvDisp_m.txt";
	std::string file3 = "InSurfEvDisp_cr.txt";
	std::string file4 = "InSurfEvDisp_cl.txt";
	std::string file5 = "InSurfEvDisp_ct.txt";
#endif

	//std::string file2 = "InSurfEvTraction.txt";
	//std::string file2 = "ouSurfEvDisp.txt";
	std::string file6 = "InSurfEvTraction.txt";
	std::string file7 = "InSurfNormailStress_m.txt";
	std::string file8 = "InSurfNormailStress_cr.txt";
	std::string file9 = "InSurfNormailStress_cl.txt";
	std::string file10 = "InSurfNormailStress_ct.txt";
	//std::ofstream of_traction(file1);
	//std::ofstream of_inSurfEvDisp(file2);
	std::ofstream of;
	//std::streampos pos1, pos2, pos7;
	if (numPlotPoints > 0)
	{
		of.open(file1);
		of << "0\n";
		of.close();

		gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, 0, numPlotPoints, meshPlot);
		//gsWriteParaviewMultiPhysics(fields, filenameParaview + "mesh", numPlotPoints, 1, 0);

#if !DISPCONTROL
		of.open(file2);
		of << "0\n";
		of.close();

		of.open(file3);
		of << "0\n";
		of.close();

		of.open(file4);
		of << "0\n";
		of.close();

		of.open(file5);
		of << "0\n";
		of.close();
#endif
		of.open(file6);
		of << "0\n";
		of.close();

		of.open(file7);
		of << "0\n";
		of.close();

		of.open(file8);
		of << "0\n";
		of.close();

		of.open(file9);
		of << "0\n";
		of.close();

		of.open(file10);
		of << "0\n";
		of.close();
	}
#if 0
	gsInfo << "Solving...\n";
	index_t numStepsPerFrame = numLoadSteps / numFrames;
	index_t cs = 0;
	index_t frame = 0;
	gsStopwatch clock;
	clock.restart();
	//bool bisectionProcess = false;
	real_t tracValue = 0.0;
	real_t SurfaceTensionValue = 0.0;
	for (int i = 0; i < numLoadSteps; i++)
	{
		gsInfo << "Step " << i + 1 << " of " << numLoadSteps;
#if DISPCONTROL
#if FIX_INNER_SURFACE
		real_t totalDisp = 0;
		assembler.options().setReal("SurfaceTension", (i + 1) * surfaceTension / numLoadSteps);
#else
		real_t totalDisp = dispValue * (i + 1);
		gsInfo << " with displacement: " << totalDisp << "\n";
#endif
		of.open(file1, std::ios_base::app);
		of << totalDisp << "\n";
		of.close();
		//dispX = gsFunctionExpr<>(util::to_string(dispValue) + "*x/" + util::to_string(ir), 3);
		//dispY = gsFunctionExpr<>(util::to_string(dispValue) + "*y/" + util::to_string(ir), 3);
		//dispZ = gsFunctionExpr<>(util::to_string(dispValue) + "*z/" + util::to_string(ir), 3);
#else
		real_t SurfaceTensionValue_old = SurfaceTensionValue;
		real_t tracValue_old = tracValue;
		if (i < numPreStep)
		{
			SurfaceTensionValue += DeltaST;
			//assembler.options().setReal("SurfaceTension", (i + 1) * surfaceTension / numPreStep);
			//gsInfo << " with SurfaceTension: " << (i + 1) * surfaceTension / numPreStep << "\n";
			assembler.options().setReal("SurfaceTension", SurfaceTensionValue);
			gsInfo << " with SurfaceTension: " << SurfaceTensionValue << "\n";

		}
		else
		{
			tracValue += DeltaTrac;
			//real_t tracValue = maxLoad / numLoadSteps * (i + 1 - numPreStep);
			gsInfo << " with traction: " << tracValue << "\n";

			traction = gsFunctionExpr<>(
				util::to_string(tracValue) + "*x/" + util::to_string(ir),
				util::to_string(tracValue) + "*y/" + util::to_string(ir),
				util::to_string(tracValue) + "*z/" + util::to_string(ir), 3);
		}
#endif
		gsMatrix<real_t> solution_old = solver.solution();
		assembler.refresh();
		solver.solve(numLoadSteps);

#if 1
		gsSparseSolver<>::SimplicialLDLT LDLTsolver(assembler.matrix());
		gsVector<> stabilityVec = LDLTsolver.vectorD();
		real_t stability = stabilityVec.colwise().minCoeff()[0];
		gsInfo << "Stability = " << stability << "\n";
		if (stability < -bisecTol)
		{
			gsInfo << "Bifurcation point found! Start to calculate bifurcation point.\n\n";
			//bisectionProcess = true;
			real_t DeltaST_bisec = DeltaST;
			real_t SurfaceTensionValue_bisec = SurfaceTensionValue_old;
			real_t DeltaTrac_bisec = DeltaTrac;
			real_t tracValue_bisec = tracValue_old;
			solver.setStepNum(1);
			index_t numBisection = 0;
			do
			{
				gsInfo << "Bisection Step " << solver.getStepNum() << " of " << maxBisecNum;
				solver.setSolutionVector(solution_old);

				if (i < numPreStep)
					DeltaST_bisec = DeltaST_bisec / 2.;
				else
					DeltaTrac_bisec = DeltaTrac_bisec / 2.;

				if (stability > 0)
				{
					//solution_old = solver.solution();
					if (i < numPreStep)
					{
						SurfaceTensionValue_old = SurfaceTensionValue_bisec;
					}
					else
					{
						tracValue_old = tracValue_bisec;
					}
				}

				if (i < numPreStep)
				{
					SurfaceTensionValue_bisec = SurfaceTensionValue_old + DeltaST_bisec;
					assembler.options().setReal("SurfaceTension", SurfaceTensionValue_bisec);
					gsInfo << " with SurfaceTension: " << SurfaceTensionValue_bisec << "\n";
				}
				else
				{
					tracValue_bisec = tracValue_old + DeltaTrac_bisec;
					traction = gsFunctionExpr<>(
						util::to_string(tracValue_bisec) + "*x/" + util::to_string(ir),
						util::to_string(tracValue_bisec) + "*y/" + util::to_string(ir),
						util::to_string(tracValue_bisec) + "*z/" + util::to_string(ir), 3);
					gsInfo << " with traction: " << tracValue_bisec << "\n";
				}

				solver.solve(maxBisecNum);
				if (solver.getStatus() == solver_status::interrupted)
					break;
				assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
				assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
				assembler.constructCauchyStressesExtension(displacement, reaction, boundary::front, stress_components::normal_3D_vector);

				if (i < numPreStep)
				{
					of.open("SurfaceTension_bisec.txt");
					of << SurfaceTensionValue_bisec << "\n";
					of.close();
				}
				else
				{
					of.open("traction_bisec.txt");
					of << tracValue_bisec << "\n";
					of.close();
				}
				gsMatrix<> p(3, 1);
				gsMatrix<> A(3, 1);
				gsMatrix<> B(3, 1);
				p << 1., 1., 0.;
				A = displacement.patch(0).eval(p);
				B = reaction.function(0).eval(p);

				of.open("InSurfEvDisp_m_bisec.txt");

				if (A.at(0) > 0)
					of << math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2)) << "\n";
				else
					of << -1 * (math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2))) << "\n";

				of.close();

				of.open("InSurfNormailStress_m_bisec.txt");
				of << math::sqrt(math::pow(B.at(0), 2) + math::pow(B.at(1), 2) + math::pow(B.at(2), 2)) << "\n";
				of.close();

				gsSparseSolver<>::SimplicialLDLT LDLTsolver(assembler.matrix());
				gsVector<> stabilityVec = LDLTsolver.vectorD();
				stability = stabilityVec.colwise().minCoeff()[0];
				gsInfo << "Stability = " << stability << "\n";
				if (abs(stability) > bisecTol)
					gsInfo << "Tolerance not reached! Bisection process continue.\n";
				numBisection++;
			} while (abs(stability) > bisecTol && numBisection < maxBisecNum);
			if (solver.getStatus() == solver_status::converged)
				gsInfo << "Bisection process finished! The critical pressure is: " << tracValue_bisec << " \n";
			else
				gsInfo << "Bisection process was interrupted!" << " \n";

			solver.setStepNum(i);

			gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, frame, numPlotPoints, meshPlot);

			break;
		}
#endif
		if (i < numPreStep)
		{
			of.open(file1, std::ios_base::app);
			of << "0\n";
			of.close();
		}
		else
		{
			of.open(file1, std::ios_base::app);
			of << tracValue << "\n";
			of.close();
		}
		assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
		assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
		assembler.constructCauchyStressesExtension(displacement, reaction, boundary::front, stress_components::normal_3D_vector);
		cs++;
		//assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
		if (numPlotPoints > 0 && cs == numStepsPerFrame)
		{
			frame++;
			gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, frame, numPlotPoints, meshPlot);
			//gsWriteHistoryOutputBoundaryResults(fields, file2, "Displacement", boundary::front, 100);
			//gsWriteHistoryOutputBoundaryResults(fields2, file7, "Reaction", boundary::front, 100);
			//gsWriteHistoryOutputBoundaryResults(fields, file3, "Displacement", boundary::back, 100);

#if !DISPCONTROL
			gsMatrix<> p(3, 1);
			gsMatrix<> A(3, 1);
			gsMatrix<> B(3, 1);

			p << 1., 1., 0.;
			A = displacement.patch(0).eval(p);
			B = reaction.function(0).eval(p);
			of.open(file2, std::ios_base::app);
			if (A.at(0) > 0)
				of << math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2)) << "\n";
			else
				of << -1 * (math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2))) << "\n";
			of.close();
			of.open(file7, std::ios_base::app);
			of << math::sqrt(math::pow(B.at(0), 2) + math::pow(B.at(1), 2) + math::pow(B.at(2), 2)) << "\n";
			of.close();

			p << 0., 0., 0.;
			A = displacement.patch(0).eval(p);
			B = reaction.function(0).eval(p);
			of.open(file3, std::ios_base::app);
			if (A.at(0) > 0)
				of << math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2)) << "\n";
			else
				of << -1 * (math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2))) << "\n";
			of.close();
			of.open(file8, std::ios_base::app);
			of << math::sqrt(math::pow(B.at(0), 2) + math::pow(B.at(1), 2) + math::pow(B.at(2), 2)) << "\n";
			of.close();

			p << 1., 0., 0.;
			A = displacement.patch(0).eval(p);
			B = reaction.function(0).eval(p);
			of.open(file4, std::ios_base::app);
			if (A.at(0) > 0)
				of << math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2)) << "\n";
			else
				of << -1 * (math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2))) << "\n";
			of.close();
			of.open(file9, std::ios_base::app);
			of << math::sqrt(math::pow(B.at(0), 2) + math::pow(B.at(1), 2) + math::pow(B.at(2), 2)) << "\n";
			of.close();

			p << 0., 1., 0.;
			A = displacement.patch(0).eval(p);
			B = reaction.function(0).eval(p);
			of.open(file5, std::ios_base::app);
			if (A.at(0) > 0)
				of << math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2)) << "\n";
			else
				of << -1 * (math::sqrt(math::pow(A.at(0), 2) + math::pow(A.at(1), 2) + math::pow(A.at(2), 2))) << "\n";
			of.close();
			of.open(file10, std::ios_base::app);
			of << math::sqrt(math::pow(B.at(0), 2) + math::pow(B.at(1), 2) + math::pow(B.at(2), 2)) << "\n";
			of.close();
#endif

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
