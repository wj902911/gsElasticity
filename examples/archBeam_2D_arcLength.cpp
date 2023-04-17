/// This is the "Cook's membrane" benchmark solved using the nonlinear elasticity solver.
/// The problem description and reference solutions can be found in the Ph.D. thesis of O.Weeger
/// "Isogeometric Finite Element Analysis of Nonlinear Structural Vibrations", 2015.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsStructuralAnalysis/gsALMBase.h>
#include <gsStructuralAnalysis/gsALMCrisfield.h>

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

    //std::string filename = ELAST_DATA_DIR"/square.xml";
    index_t material = 2;
    real_t youngsModulus = 10.0;
    real_t poissonsRatio = 0.4;
    real_t surfaceTension = 1.0;
    real_t surfaceYoungsModulus = 0.0;
    real_t surfacePoissonsRatio = 0.4;
    index_t numUniRef_l = 2;
    index_t numUniRef_t = 1;
    index_t numDegElev = 1;
    index_t numPlotPoints = 1000;
    index_t numLoadSteps = 100;
    //index_t numPreSteps = 10;
    index_t maxIter = 20;
    real_t maxLoad = -1.0;

    real_t dL = 0.05; // General arc length
    real_t dLb = 0.05; // Ard length to find bifurcation
    bool adaptive = true;
    real_t tau = 1e4;
    real_t tol = 1e-6;
    real_t tolU = 1e-6;
    real_t tolF = 1e-3;
    real_t relax = 1.0;

    bool quasiNewton = false;
    int quasiNewtonInt = -1;
    bool SingularPoint = true;


    // minimalistic user interface for terminal
    gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
    cmd.addReal("p", "poisson", "Poisson's ratio used in the material law", poissonsRatio);
    cmd.addReal("v","sufacePoisson","Poisson's ratio used in the surface material law", surfacePoissonsRatio);
    cmd.addInt("r","refine","Number of uniform refinement application", numUniRef_l);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addInt("s","point","Number of points to plot to Paraview",numPlotPoints);
    cmd.addReal("t", "surfTen", "Constant surface tension", surfaceTension);
    cmd.addInt("n", "numSteps", "Number of load steps", numLoadSteps);
    cmd.addInt("m", "maxIter", "Max iteration number", maxIter);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    gsKnotVector<> KV1(0, 1, 0, 3);
    gsKnotVector<> KV2(0, 1, 0, 2);
	
    gsMatrix<real_t> C(6, 2);
    gsMatrix<real_t> W(6, 1);
    gsMatrix<real_t> C1(3, 2);
    gsMatrix<real_t> W1(3, 1);

    real_t s3 = math::sqrt(3.);

	C1 << -0.5,      s3 / 2.,
		    0., 2. * s3 / 3.,
		   0.5,      s3 / 2.;
	W1 << 1., s3 / 2., 1.;
	
    C << C1, C1 * 1.05;
    W << W1, W1;
    gsTensorNurbs<2, real_t> patch1(KV1, KV2, C, W);

    // scanning geometry
    gsMultiPatch<> geometry;
    geometry.addPatch(patch1);
    geometry.computeTopology();
    // creating bases
    gsBSplineBasis<real_t>* Bu = new gsBSplineBasis<real_t>(KV1);
    gsBSplineBasis<real_t>* Bv = new gsBSplineBasis<real_t>(KV2);
    for (int i = 0; i < numUniRef_l; i++)
    {
        gsSparseMatrix<real_t, RowMajor> B[2];
        gsSparseMatrix<real_t, RowMajor> trans;
        Bu->uniformRefine_withTransfer(B[0], 1);
        Bv->uniformRefine_withTransfer(B[1], 0);
        tensorCombineTransferMatrices<2, real_t>(B, trans);
        W = trans * W;
    }
    for (int i = 0; i < numUniRef_t; i++)
    {
        gsSparseMatrix<real_t, RowMajor> B[2];
        gsSparseMatrix<real_t, RowMajor> trans;
        Bu->uniformRefine_withTransfer(B[0], 0);
        Bv->uniformRefine_withTransfer(B[1], 1);
        tensorCombineTransferMatrices<2, real_t>(B, trans);
        W = trans * W;
    }
    gsTensorBSplineBasis<2, real_t>* tbasis = new gsTensorBSplineBasis<2, real_t>(Bu, Bv);
    gsTensorNurbsBasis<2, real_t>* basis = new gsTensorNurbsBasis<2, real_t>(tbasis, give(W));
    gsMultiBasis<> basisDisplacement;
    basisDisplacement.setTopology(geometry.topology());
    for (int i = 0; i < geometry.nPatches(); i++)
        basisDisplacement.addBasis(basis->clone());
	
    for (index_t i = 0; i < 1; ++i)
        basisDisplacement.degreeElevate(1, 1);
    for (index_t i = 0; i < numDegElev; ++i)
        basisDisplacement.degreeElevate();
    //for (index_t i = 0; i < numUniRef; ++i)
        //basisDisplacement.uniformRefine();
    //gsGeometry<>* pGeom = &basisDisplacement;


    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // neumann BC
    gsConstantFunction<> f(0., maxLoad, 2);
	// displacement BC
    //gsConstantFunction<> f(0.1, 2);
    //gsConstantFunction<> d(-2.7e-3, 2);
    // elasticSurface BC
    //gsConstantFunction<> surfaceTension(4., 1);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    bcInfo.addCornerValue(boundary::southwest, 0, 0, 0);
    bcInfo.addCornerValue(boundary::southwest, 0, 0, 1);
    bcInfo.addCornerValue(boundary::southeast, 0, 0, 0);
    bcInfo.addCornerValue(boundary::southeast, 0, 0, 1);
    //for (index_t d = 0; d < 2; ++d)
    //{
        //bcInfo.addCornerValue(boundary::southwest, 0, 0, d);
        //bcInfo.addCornerValue(boundary::northwest, 0, 0, d);
        //bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
    //}
    bcInfo.addCondition(0, boundary::north, condition_type::neumann, &f);
    //bcInfo.addCondition(0, boundary::east, condition_type::dirichlet, &d, 1);
    //bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, 0);
    //bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, 1);
    //bcInfo.addCondition(0, boundary::south, condition_type::robin, nullptr);
    //bcInfo.addCondition(0, boundary::east, condition_type::robin, nullptr);

    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler
    gsElasticityAssemblerElasticSurface<real_t> assembler(geometry, basisDisplacement, bcInfo, g);
    assembler.options().setReal("YoungsModulus", youngsModulus);
    assembler.options().setReal("PoissonsRatio", poissonsRatio);
    assembler.options().setReal("SurfaceYoungsModulus", surfaceYoungsModulus);
    assembler.options().setReal("SurfacePoissonsRatio", surfacePoissonsRatio);
    assembler.options().setReal("SurfaceTension", surfaceTension);
    assembler.options().setInt("MaterialLaw", material_law::hooke);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    gsStopwatch stopwatch;
    real_t time = 0.0;
	
    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
    typedef std::function<gsSparseMatrix<real_t>(gsVector<real_t> const&)>                                Jacobian_t;
    typedef std::function<gsVector<real_t>(gsVector<real_t> const&, real_t, gsVector<real_t> const&) >   ALResidual_t;
    Jacobian_t Jacobian = [&time, &stopwatch, &assembler, &fixedDofs](gsVector<real_t> const& x)
    {
        stopwatch.restart();
        assembler.assemble(x, fixedDofs);
        time += stopwatch.stop();

        gsSparseMatrix<real_t> m = assembler.matrix();
        // gsInfo<<"matrix = \n"<<m.toDense()<<"\n";
        return m;
    };
    // Function for the Residual
    ALResidual_t ALResidual = [&time, &stopwatch, &assembler, &fixedDofs](gsVector<real_t> const& x, real_t lam, gsVector<real_t> const& force)
    {
        stopwatch.restart();
        assembler.assemble(x, fixedDofs);
        gsVector<real_t> Fint = -(assembler.rhs() - force);
        gsVector<real_t> result = Fint - lam * force;
        time += stopwatch.stop();
        return result; // - lam * force;
    };
    // Assemble linear system to obtain the force vector
    assembler.assemble();
    gsVector<> Force = assembler.rhs();
    assembler.options().setInt("MaterialLaw", material_law::neo_hooke_ln);
    gsALMBase<real_t>* arcLength;
    arcLength = new gsALMCrisfield<real_t>(Jacobian, ALResidual, Force);
	
    arcLength->options().setString("Solver", "SimplicialLDLT"); // LDLT solver
    arcLength->options().setInt("BifurcationMethod", 0); // 0: determinant, 1: eigenvalue
    arcLength->options().setReal("Length", dLb);
    arcLength->options().setInt("AngleMethod", 0); // 0: step, 1: iteration
    arcLength->options().setSwitch("AdaptiveLength", adaptive);
    arcLength->options().setInt("AdaptiveIterations", 6);
    arcLength->options().setReal("MaxLength", 0.05);
    arcLength->options().setReal("Perturbation", tau);
    arcLength->options().setReal("Scaling", 1.0);
    arcLength->options().setReal("Tol", tol);
    arcLength->options().setReal("TolU", tolU);
    arcLength->options().setReal("TolF", tolF);
    arcLength->options().setInt("MaxIter", maxIter);
    arcLength->options().setSwitch("Verbose", true);
    arcLength->options().setReal("Relaxation", relax);
    if (quasiNewtonInt > 0)
    {
        quasiNewton = true;
        arcLength->options().setInt("QuasiIterations", quasiNewtonInt);
    }
    arcLength->options().setSwitch("Quasi", quasiNewton);
    arcLength->applyOptions();
    arcLength->initialize();

    gsMatrix<> solVector = arcLength->solutionU();
    gsMultiPatch<> displacement;
    assembler.constructSolution(solVector, fixedDofs, displacement);
    gsPiecewiseFunction<> stresses;
    assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
    // constructing an IGA field (geometry + solution)
    gsField<> displacementField(assembler.patches(), displacement);
    //gsField<> meshAndCnet(displacement, displacement);
    gsField<> stressField(assembler.patches(), stresses, true);
    std::map<std::string, const gsField<>*> fields;
    fields["Displacement"] = &displacementField;
    fields["von Mises"] = &stressField;
    std::string filenameParaview = "ArchBeam_ALM_" + util::to_string(numUniRef_l) + "_" + util::to_string(numUniRef_t) + "_" + util::to_string(youngsModulus) + "_" + util::to_string(surfaceYoungsModulus) + "_";
    gsParaviewCollection collection(filenameParaview);
    std::string file1 = "traction.txt";
    std::string file2 = "displacement.txt";
    std::ofstream of;

    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, 0, numPlotPoints, true);

        of.open(file1);
        of << "0\n";
        of.close();

        of.open(file2);
        of << "0\n";
        of.close();
    }
	
#if 1
    // Make objects for previous solutions
    real_t Lold = 0;
    gsMatrix<> Uold = Force;
    Uold.setZero();
	
    real_t indicator = 0.0;
    arcLength->setIndicator(indicator); // RESET INDICATOR
    //bool bisected = false;
    real_t dLb0 = dLb;
    for (index_t k = 0; k < numLoadSteps; k++)
    {
        gsInfo << "Load step " << k << "\n";
        arcLength->step();
		
        if (!(arcLength->converged()))
        {
            gsInfo << "Error: Loop terminated, arc length method did not converge.\n";
            //dLb = dLb / 2.;
            real_t currentLength = arcLength->getLength();
            arcLength->setLength(currentLength / 2.);
            arcLength->setSolution(Uold, Lold);
            //bisected = true;
            k -= 1;
            continue;
        }

        if (SingularPoint)
        {
            arcLength->computeStability(arcLength->solutionU(), quasiNewton);
            if (arcLength->stabilityChange())
            {
                gsInfo << "Bifurcation spotted!" << "\n";
                arcLength->computeSingularPoint(1e-4, 5, Uold, Lold, 1e-10, 1e-1, false);
                arcLength->switchBranch();
                dLb0 = dLb = dL;
                arcLength->setLength(dLb);
            }
        }
		
        indicator = arcLength->indicator();
        solVector = arcLength->solutionU();
        Uold = solVector;
        //gsInfo << Uold << "\n\n";
        Lold = arcLength->solutionL();
        //gsInfo << Lold << "\n\n";
        assembler.constructSolution(solVector, fixedDofs, displacement);
        assembler.constructCauchyStresses(displacement, stresses, stress_components::von_mises);
		
        gsInfo << "Total ellapsed assembly time: " << time << " s\n";

        if (numPlotPoints > 0)
        {
            gsWriteParaviewMultiPhysicsTimeStepWithMesh(fields, filenameParaview, collection, k + 1, numPlotPoints, true);

            of.open(file1, std::ios_base::app);
            of << -Lold * maxLoad << "\n";
            of.close();
			
            gsMatrix<> A(2, 1);
            A << .5, .5;
            A = displacement.patch(0).eval(A);

            of.open(file2, std::ios_base::app);
            of << -1 * A.at(1) << "\n";
            of.close();
        }
		
#if 0
        if (!bisected)
        {
            dLb = dLb0;
            arcLength->setLength(dLb);
        }
        bisected = false;
#endif	
    }
#endif

    delete arcLength;

    if (numPlotPoints > 0)
    {
        collection.save();
    }

    return 0;
}
