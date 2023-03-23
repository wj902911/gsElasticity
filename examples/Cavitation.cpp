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
#include <gsNurbs/gsNurbsCreator.h>

using namespace gismo;

gsMatrix<real_t> generateRotationMatrix(real_t angle, gsVector<index_t> rotationAxis)
{
	gsMatrix<real_t> rotationMatrix(3, 3);
	rotationMatrix << cos(angle) + rotationAxis[0] * rotationAxis[0] * (1 - cos(angle)), rotationAxis[0] * rotationAxis[1] * (1 - cos(angle)) - rotationAxis[2] * sin(angle), rotationAxis[0] * rotationAxis[2] * (1 - cos(angle)) + rotationAxis[1] * sin(angle),
		rotationAxis[1] * rotationAxis[0] * (1 - cos(angle)) + rotationAxis[2] * sin(angle), cos(angle) + rotationAxis[1] * rotationAxis[1] * (1 - cos(angle)), rotationAxis[1] * rotationAxis[2] * (1 - cos(angle)) - rotationAxis[0] * sin(angle),
		rotationAxis[2] * rotationAxis[0] * (1 - cos(angle)) - rotationAxis[1] * sin(angle), rotationAxis[2] * rotationAxis[1] * (1 - cos(angle)) + rotationAxis[0] * sin(angle), cos(angle) + rotationAxis[2] * rotationAxis[2] * (1 - cos(angle));
	return rotationMatrix;
}

int main(int argc, char* argv[])
{

	//=====================================//
				// Input //
	//=====================================//
	
	//geometry
	real_t r = 1.0;
	real_t ir = 0.1;
	real_t x = 0.0;
	real_t y = 0.0;
	real_t z = 0.0;

	//material
	real_t youngsModulus = 10.0;
	real_t poissonsRatio = 0.4;
	real_t surfaceTension = 0.0;
	real_t surfaceYoungsModulus = 0.0;
	real_t surfacePoissonsRatio = 0.4;

	//mesh
	index_t numUniRef = 0;
	index_t numDegElev = 0;
	index_t numDegRedu = 1;

	//solver
	index_t numLoadSteps = 50;
	index_t maxIter = 100;

	//plot
	index_t numFrames = 50;
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

	//=============================================//
		// Createing geometry and bases //
	//=============================================//

	gsKnotVector<> KV1(0, 1, 0, 5);
	gsKnotVector<> KV2(0, 1, 0, 5);
	gsKnotVector<> KV3(0, 1, 0, 2);
	gsMatrix<real_t> C(50, 3);
	gsMatrix<real_t> C1(25, 3);
	gsMatrix<real_t> C2(25, 3);
	gsMatrix<real_t> CTemp(25, 3);
	gsMatrix<real_t> W(50, 1);

	gsVector<index_t> rotAxis(3);
	rotAxis << 0, 0, 1;
	gsMatrix<real_t> rotMat = generateRotationMatrix(M_PI / 4., rotAxis);

	W <<
		5.0717967697244912, 4.5200421036033447, 4.3572655899081640, 4.5200421036033447, 5.0717967697244912,
		4.5200421036033447, 3.8660254037844384, 3.6449237056739161, 3.8660254037844384, 4.5200421036033447,
		4.3572655899081640, 3.6449237056739161, 3.4045573501530604, 3.6449237056739161, 4.3572655899081640,
		4.5200421036033447, 3.8660254037844384, 3.6449237056739161, 3.8660254037844384, 4.5200421036033447,
		5.0717967697244912, 4.5200421036033447, 4.3572655899081640, 4.5200421036033447, 5.0717967697244912,
		5.0717967697244912, 4.5200421036033447, 4.3572655899081640, 4.5200421036033447, 5.0717967697244912,
		4.5200421036033447, 3.8660254037844384, 3.6449237056739161, 3.8660254037844384, 4.5200421036033447,
		4.3572655899081640, 3.6449237056739161, 3.4045573501530604, 3.6449237056739161, 4.3572655899081640,
		4.5200421036033447, 3.8660254037844384, 3.6449237056739161, 3.8660254037844384, 4.5200421036033447,
		5.0717967697244912, 4.5200421036033447, 4.3572655899081640, 4.5200421036033447, 5.0717967697244912;

	//bottom
	CTemp <<
		-0.57735026918962562, -0.57735026918962562, -0.57735026918962562, 
		-0.31287619229159269, -0.70958730763864231, -0.70958730763864231, 
		 0.00000000000000000, -0.75402078491454050, -0.75402078491454050,
		 0.31287619229159269, -0.70958730763864231, -0.70958730763864231, 
		 0.57735026918962562, -0.57735026918962562, -0.57735026918962562, 
		
		-0.70958730763864231, -0.31287619229159269, -0.70958730763864231,
		-0.41336412579931964, -0.41336412579931964, -1.00000000000000000, 
		 0.00000000000000000, -0.45730408049064641, -1.12004618869897940, 
		 0.41336412579931964, -0.41336412579931964, -1.00000000000000000,
		 0.70958730763864231, -0.31287619229159269, -0.70958730763864231,
		
		-0.75402078491454050,  0.00000000000000000, -0.75402078491454050,
		-0.45730408049064641,  0.00000000000000000, -1.12004618869897940,
		 0.00000000000000000,  0.00000000000000000, -1.27983321817512400,
		 0.45730408049064641,  0.00000000000000000, -1.12004618869897940,
		 0.75402078491454050,  0.00000000000000000, -0.75402078491454050,
		
		-0.70958730763864231,  0.31287619229159269, -0.70958730763864231,
		-0.41336412579931964,  0.41336412579931964, -1.00000000000000000,
		 0.00000000000000000,  0.45730408049064641, -1.12004618869897940,
		 0.41336412579931964,  0.41336412579931964, -1.00000000000000000,
		 0.70958730763864231,  0.31287619229159269, -0.70958730763864231,
		
		-0.57735026918962562,  0.57735026918962562, -0.57735026918962562,
		-0.31287619229159269,  0.70958730763864231, -0.70958730763864231,
		 0.00000000000000000,  0.75402078491454050, -0.75402078491454050,
		 0.31287619229159269,  0.70958730763864231, -0.70958730763864231,
		 0.57735026918962562,  0.57735026918962562, -0.57735026918962562;
	
	C1 << CTemp * r * rotMat.transpose();
	C2 << CTemp * ir * rotMat.transpose();
	C << C2, C1;
	
	C.col(0).array() += x;
	C.col(1).array() += y;
	C.col(2).array() += z;

	gsTensorNurbs<3, real_t> patch1(KV1, KV2, KV3, C, W);
	
	//top
	CTemp <<
		-0.57735026918962562, -0.57735026918962562,  0.57735026918962562, 
		-0.31287619229159269, -0.70958730763864231,  0.70958730763864231, 
		 0.00000000000000000, -0.75402078491454050,  0.75402078491454050,
		 0.31287619229159269, -0.70958730763864231,  0.70958730763864231, 
		 0.57735026918962562, -0.57735026918962562,  0.57735026918962562, 
		
		-0.70958730763864231, -0.31287619229159269,  0.70958730763864231,
		-0.41336412579931964, -0.41336412579931964,  1.00000000000000000, 
		 0.00000000000000000, -0.45730408049064641,  1.12004618869897940, 
		 0.41336412579931964, -0.41336412579931964,  1.00000000000000000,
		 0.70958730763864231, -0.31287619229159269,  0.70958730763864231,
		
		-0.75402078491454050,  0.00000000000000000,  0.75402078491454050,
		-0.45730408049064641,  0.00000000000000000,  1.12004618869897940,
		 0.00000000000000000,  0.00000000000000000,  1.27983321817512400,
		 0.45730408049064641,  0.00000000000000000,  1.12004618869897940,
		 0.75402078491454050,  0.00000000000000000,  0.75402078491454050,
		
		-0.70958730763864231,  0.31287619229159269,  0.70958730763864231,
		-0.41336412579931964,  0.41336412579931964,  1.00000000000000000,
		 0.00000000000000000,  0.45730408049064641,  1.12004618869897940,
		 0.41336412579931964,  0.41336412579931964,  1.00000000000000000,
		 0.70958730763864231,  0.31287619229159269,  0.70958730763864231,
		
		-0.57735026918962562,  0.57735026918962562,  0.57735026918962562,
		-0.31287619229159269,  0.70958730763864231,  0.70958730763864231,
		 0.00000000000000000,  0.75402078491454050,  0.75402078491454050,
		 0.31287619229159269,  0.70958730763864231,  0.70958730763864231,
		 0.57735026918962562,  0.57735026918962562,  0.57735026918962562;

	
	C1 << CTemp * r * rotMat.transpose();
	C2 << CTemp * ir * rotMat.transpose();
	C << C2, C1;

	C.col(0).array() += x;
	C.col(1).array() += y;
	C.col(2).array() += z;

	gsTensorNurbs<3, real_t> patch2(KV1, KV2, KV3, C, W);
	
	//front u=0 yz:ux=0 || u=1 xz:uy=0
	CTemp <<
		-0.57735026918962562, -0.57735026918962562, -0.57735026918962562, 
		-0.31287619229159269, -0.70958730763864231, -0.70958730763864231, 
		 0.00000000000000000, -0.75402078491454050, -0.75402078491454050,
		 0.31287619229159269, -0.70958730763864231, -0.70958730763864231, 
		 0.57735026918962562, -0.57735026918962562, -0.57735026918962562, 
		
		-0.70958730763864231, -0.70958730763864231, -0.31287619229159269,
		-0.41336412579931964, -1.00000000000000000, -0.41336412579931964,
		 0.00000000000000000, -1.12004618869897940, -0.45730408049064641,
		 0.41336412579931964, -1.00000000000000000, -0.41336412579931964,
		 0.70958730763864231, -0.70958730763864231, -0.31287619229159269,
		
		-0.75402078491454050, -0.75402078491454050,  0.00000000000000000,
		-0.45730408049064641, -1.12004618869897940,  0.00000000000000000,
		 0.00000000000000000, -1.27983321817512400,  0.00000000000000000,
		 0.45730408049064641, -1.12004618869897940,  0.00000000000000000,
		 0.75402078491454050, -0.75402078491454050,  0.00000000000000000,
		
		-0.70958730763864231, -0.70958730763864231,  0.31287619229159269,
		-0.41336412579931964, -1.00000000000000000,  0.41336412579931964,
		 0.00000000000000000, -1.12004618869897940,  0.45730408049064641,
		 0.41336412579931964, -1.00000000000000000,  0.41336412579931964,
		 0.70958730763864231, -0.70958730763864231,  0.31287619229159269,
		
		-0.57735026918962562, -0.57735026918962562,  0.57735026918962562,
		-0.31287619229159269, -0.70958730763864231,  0.70958730763864231,
		 0.00000000000000000, -0.75402078491454050,  0.75402078491454050,
		 0.31287619229159269, -0.70958730763864231,  0.70958730763864231,
		 0.57735026918962562, -0.57735026918962562,  0.57735026918962562;

	
	C1 << CTemp * r * rotMat.transpose();
	C2 << CTemp * ir * rotMat.transpose();
	C << C2, C1;

	C.col(0).array() += x;
	C.col(1).array() += y;
	C.col(2).array() += z;

	gsTensorNurbs<3, real_t> patch3(KV1, KV2, KV3, C, W);
	
	//back u=0 xz:uy=0 || u=1 yz:ux=0
	CTemp <<
		-0.57735026918962562,  0.57735026918962562, -0.57735026918962562, 
		-0.31287619229159269,  0.70958730763864231, -0.70958730763864231, 
		 0.00000000000000000,  0.75402078491454050, -0.75402078491454050,
		 0.31287619229159269,  0.70958730763864231, -0.70958730763864231, 
		 0.57735026918962562,  0.57735026918962562, -0.57735026918962562, 
		
		-0.70958730763864231,  0.70958730763864231, -0.31287619229159269,
		-0.41336412579931964,  1.00000000000000000, -0.41336412579931964,
		 0.00000000000000000,  1.12004618869897940, -0.45730408049064641,
		 0.41336412579931964,  1.00000000000000000, -0.41336412579931964,
		 0.70958730763864231,  0.70958730763864231, -0.31287619229159269,
		
		-0.75402078491454050,  0.75402078491454050,  0.00000000000000000,
		-0.45730408049064641,  1.12004618869897940,  0.00000000000000000,
		 0.00000000000000000,  1.27983321817512400,  0.00000000000000000,
		 0.45730408049064641,  1.12004618869897940,  0.00000000000000000,
		 0.75402078491454050,  0.75402078491454050,  0.00000000000000000,
		
		-0.70958730763864231,  0.70958730763864231,  0.31287619229159269,
		-0.41336412579931964,  1.00000000000000000,  0.41336412579931964,
		 0.00000000000000000,  1.12004618869897940,  0.45730408049064641,
		 0.41336412579931964,  1.00000000000000000,  0.41336412579931964,
		 0.70958730763864231,  0.70958730763864231,  0.31287619229159269,
		
		-0.57735026918962562,  0.57735026918962562,  0.57735026918962562,
		-0.31287619229159269,  0.70958730763864231,  0.70958730763864231,
		 0.00000000000000000,  0.75402078491454050,  0.75402078491454050,
		 0.31287619229159269,  0.70958730763864231,  0.70958730763864231,
		 0.57735026918962562,  0.57735026918962562,  0.57735026918962562;

	
	C1 << CTemp * r * rotMat.transpose();
	C2 << CTemp * ir * rotMat.transpose();
	C << C2, C1;

	C.col(0).array() += x;
	C.col(1).array() += y;
	C.col(2).array() += z;

	gsTensorNurbs<3, real_t> patch4(KV1, KV2, KV3, C, W);
	//left v=0 yz:ux=0 || v=1 xz:uy=0
	CTemp <<
		-0.57735026918962562, -0.57735026918962562,  0.57735026918962562,
		-0.70958730763864231, -0.70958730763864231,  0.31287619229159269,
		-0.75402078491454050, -0.75402078491454050,  0.00000000000000000,
		-0.70958730763864231, -0.70958730763864231, -0.31287619229159269,
		-0.57735026918962562, -0.57735026918962562, -0.57735026918962562,
		
		-0.70958730763864231, -0.31287619229159269,  0.70958730763864231,
		-1.00000000000000000, -0.41336412579931964,  0.41336412579931964,
		-1.12004618869897940, -0.45730408049064641,  0.00000000000000000,
		-1.00000000000000000, -0.41336412579931964, -0.41336412579931964,
		-0.70958730763864231, -0.31287619229159269, -0.70958730763864231,
		
		-0.75402078491454050,  0.00000000000000000,  0.75402078491454050,
		-1.12004618869897940,  0.00000000000000000,  0.45730408049064641,
		-1.27983321817512400,  0.00000000000000000,  0.00000000000000000,
		-1.12004618869897940,  0.00000000000000000, -0.45730408049064641,
		-0.75402078491454050,  0.00000000000000000, -0.75402078491454050,
		
		-0.70958730763864231,  0.31287619229159269,  0.70958730763864231,
		-1.00000000000000000,  0.41336412579931964,  0.41336412579931964,
		-1.12004618869897940,  0.45730408049064641,  0.00000000000000000,
		-1.00000000000000000,  0.41336412579931964, -0.41336412579931964,
		-0.70958730763864231,  0.31287619229159269, -0.70958730763864231,
		
		-0.57735026918962562,  0.57735026918962562,  0.57735026918962562,
		-0.70958730763864231,  0.70958730763864231,  0.31287619229159269,
		-0.75402078491454050,  0.75402078491454050,  0.00000000000000000,
		-0.70958730763864231,  0.70958730763864231, -0.31287619229159269,
		-0.57735026918962562,  0.57735026918962562, -0.57735026918962562;

	
	C1 << CTemp * r * rotMat.transpose();
	C2 << CTemp * ir * rotMat.transpose();
	C << C2, C1;

	C.col(0).array() += x;
	C.col(1).array() += y;
	C.col(2).array() += z;

	gsTensorNurbs<3, real_t> patch5(KV1, KV2, KV3, C, W);

	//right v=0 xz:uy=0 || v=1 yz:ux=0
	CTemp <<
		0.57735026918962562, -0.57735026918962562,  0.57735026918962562,
		0.70958730763864231, -0.70958730763864231,  0.31287619229159269,
		0.75402078491454050, -0.75402078491454050,  0.00000000000000000,
		0.70958730763864231, -0.70958730763864231, -0.31287619229159269,
		0.57735026918962562, -0.57735026918962562, -0.57735026918962562,
		
		0.70958730763864231, -0.31287619229159269,  0.70958730763864231,
		1.00000000000000000, -0.41336412579931964,  0.41336412579931964,
		1.12004618869897940, -0.45730408049064641,  0.00000000000000000,
		1.00000000000000000, -0.41336412579931964, -0.41336412579931964,
		0.70958730763864231, -0.31287619229159269, -0.70958730763864231,
		
		0.75402078491454050,  0.00000000000000000,  0.75402078491454050,
		1.12004618869897940,  0.00000000000000000,  0.45730408049064641,
		1.27983321817512400,  0.00000000000000000,  0.00000000000000000,
		1.12004618869897940,  0.00000000000000000, -0.45730408049064641,
		0.75402078491454050,  0.00000000000000000, -0.75402078491454050,
		
		0.70958730763864231,  0.31287619229159269,  0.70958730763864231,
		1.00000000000000000,  0.41336412579931964,  0.41336412579931964,
		1.12004618869897940,  0.45730408049064641,  0.00000000000000000,
		1.00000000000000000,  0.41336412579931964, -0.41336412579931964,
		0.70958730763864231,  0.31287619229159269, -0.70958730763864231,
		
		0.57735026918962562,  0.57735026918962562,  0.57735026918962562,
		0.70958730763864231,  0.70958730763864231,  0.31287619229159269,
		0.75402078491454050,  0.75402078491454050,  0.00000000000000000,
		0.70958730763864231,  0.70958730763864231, -0.31287619229159269,
		0.57735026918962562,  0.57735026918962562, -0.57735026918962562;

	
	C1 << CTemp * r * rotMat.transpose();
	C2 << CTemp * ir * rotMat.transpose();
	C << C2, C1;

	C.col(0).array() += x;
	C.col(1).array() += y;
	C.col(2).array() += z;

	gsTensorNurbs<3, real_t> patch6(KV1, KV2, KV3, C, W);
	
	gsMultiPatch<> geometry;
	geometry.addPatch(patch1);//-z
	geometry.addPatch(patch2);//+z
	geometry.addPatch(patch3);//+x,-y u=0 yz:ux=0 || u=1 xz:uy=0
	geometry.addPatch(patch4);//-x,+y u=0 xz:uy=0 || u=1 yz:ux=0
	geometry.addPatch(patch5);//-x,-y v=0 yz:ux=0 || v=1 xz:uy=0
	geometry.addPatch(patch6);//+x,+y v=0 xz:uy=0 || v=1 yz:ux=0

	geometry.computeTopology();
	
	//gsFileData<> fd;
	//fd << geometry;
	//std::string output("geo");
	//fd.save(output);
	//gsInfo << geometry << "\n";

	// creating bases
	gsMultiBasis<> basisDisplacement(geometry);
	for (index_t i = 0; i < numDegElev; ++i)
		basisDisplacement.degreeElevate();
	for (index_t i = 0; i < numDegRedu; ++i)
		basisDisplacement.degreeReduce();
	for (index_t i = 0; i < numUniRef; ++i)
		basisDisplacement.uniformRefine();

	//=============================================//
		// Setting loads and boundary conditions //
	//=============================================//
	
	gsBoundaryConditions<> bcInfo;
	
	bcInfo.addCornerValue(boundary::southwestback, 0, 0, 0);
	//bcInfo.addCornerValue(boundary::southwestback, 0, 0, 2);
	bcInfo.addCornerValue(boundary::southeastback, 0, 0, 1);
	//bcInfo.addCornerValue(boundary::southeastback, 0, 0, 2);
	bcInfo.addCornerValue(boundary::northwestback, 0, 0, 1);
	//bcInfo.addCornerValue(boundary::northwestback, 0, 0, 2);
	bcInfo.addCornerValue(boundary::northeastback, 0, 0, 0);
	//bcInfo.addCornerValue(boundary::northeastback, 0, 0, 2);
	//bcInfo.addCondition(0, boundary::southeastback, condition_type::dirichlet, nullptr, 2);
	//bcInfo.addCondition(0, boundary::northwestback, condition_type::dirichlet, nullptr, 2);
	//bcInfo.addCondition(0, boundary::back, condition_type::dirichlet, nullptr, 0);
	//bcInfo.addCondition(0, boundary::back, condition_type::dirichlet, nullptr, 1);
	//bcInfo.addCondition(0, boundary::back, condition_type::dirichlet, nullptr, 2);
	//bcInfo.addCondition(2, boundary::west, condition_type::dirichlet, nullptr, 0);
	//bcInfo.addCondition(2, boundary::east, condition_type::dirichlet, nullptr, 1);
	//bcInfo.addCondition(3, boundary::west, condition_type::dirichlet, nullptr, 1);
	//bcInfo.addCondition(3, boundary::east, condition_type::dirichlet, nullptr, 0);

	// neumann BC
	real_t s2 = sqrt(2.);
	gsConstantFunction<> f1(0., 0., -1., 3);
	gsConstantFunction<> f2(0., 0., 1., 3);
	gsConstantFunction<> f3(s2 / .2, -s2 / .2, 0., 3);
	gsConstantFunction<> f4(-s2 / .2, s2 / .2, 0., 3);
	gsConstantFunction<> f5(-s2 / .2, -s2 / .2, 0., 3);
	gsConstantFunction<> f6(s2 / .2, s2 / .2, 0., 3);
	bcInfo.addCondition(0, boundary::front, condition_type::neumann, &f1);
	bcInfo.addCondition(1, boundary::front, condition_type::neumann, &f2);
	//bcInfo.addCondition(2, boundary::front, condition_type::neumann, &f3);
	//bcInfo.addCondition(3, boundary::front, condition_type::neumann, &f4);
	//bcInfo.addCondition(4, boundary::front, condition_type::neumann, &f5);
	//bcInfo.addCondition(5, boundary::front, condition_type::neumann, &f6);

	// source function, rhs
	gsConstantFunction<> g(0., 0., 0., 3);

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
	std::string filenameParaview = "Cavitation_" + util::to_string(numUniRef) + "_" + util::to_string(youngsModulus) + "_" + util::to_string(surfaceYoungsModulus) + "_";
	gsParaviewCollection collection(filenameParaview);

	if (numPlotPoints > 0)
	{
		gsWriteParaviewMultiPhysicsTimeStep(fields, filenameParaview, collection, 0, numPlotPoints);
		gsWriteParaviewMultiPhysics(fields, filenameParaview + "mesh", numPlotPoints, 1, 0);
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
		const gsVector<real_t> val1 = Eigen::Vector3d(0., 0., -1. / numLoadSteps * (i + 1));
		const gsVector<real_t> val2 = Eigen::Vector3d(0., 0., 1. / numLoadSteps * (i + 1));
		const gsVector<real_t> val3 = Eigen::Vector3d(s2 / .2 / numLoadSteps * (i + 1), -s2 / .2 / numLoadSteps * (i + 1), 0.);
		const gsVector<real_t> val4 = Eigen::Vector3d(-s2 / .2 / numLoadSteps * (i + 1), s2 / .2 / numLoadSteps * (i + 1), 0.);
		const gsVector<real_t> val5 = Eigen::Vector3d(-s2 / .2 / numLoadSteps * (i + 1), -s2 / .2 / numLoadSteps * (i + 1), 0.);
		const gsVector<real_t> val6 = Eigen::Vector3d(s2 / .2 / numLoadSteps * (i + 1), s2 / .2 / numLoadSteps * (i + 1), 0.);
		f1.setValue(val1, 3);
		f2.setValue(val2, 3);
		f3.setValue(val3, 3);
		f4.setValue(val4, 3);
		f5.setValue(val5, 3);
		f6.setValue(val6, 3);
		solver.solve(numLoadSteps);
		assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
		//assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
		cs++;
		if (numPlotPoints > 0 && cs == numStepsPerFrame)
		{
			frame++;
			gsWriteParaviewMultiPhysicsTimeStep(fields, filenameParaview, collection, frame, numPlotPoints);
			cs = 0;
		}
	}
	gsInfo << "Solved the system in " << clock.stop() << "s.\n";
#endif
	if (numPlotPoints > 0)
	{
		collection.save();
	}
	return 0;
}
