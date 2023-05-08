/// This is the "Cook's membrane" benchmark solved using the nonlinear elasticity solver.
/// The problem description and reference solutions can be found in the Ph.D. thesis of O.Weeger
/// "Isogeometric Finite Element Analysis of Nonlinear Structural Vibrations", 2015.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler_elasticSurface.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

gsMatrix<real_t> generateRotationMatrix(real_t angle, gsVector<real_t> rotationAxis)
{
    gsMatrix<real_t> rotationMatrix(3, 3);
    rotationMatrix << cos(angle) + rotationAxis[0] * rotationAxis[0] * (1 - cos(angle)), rotationAxis[0] * rotationAxis[1] * (1 - cos(angle)) - rotationAxis[2] * sin(angle), rotationAxis[0] * rotationAxis[2] * (1 - cos(angle)) + rotationAxis[1] * sin(angle),
        rotationAxis[1] * rotationAxis[0] * (1 - cos(angle)) + rotationAxis[2] * sin(angle), cos(angle) + rotationAxis[1] * rotationAxis[1] * (1 - cos(angle)), rotationAxis[1] * rotationAxis[2] * (1 - cos(angle)) - rotationAxis[0] * sin(angle),
        rotationAxis[2] * rotationAxis[0] * (1 - cos(angle)) - rotationAxis[1] * sin(angle), rotationAxis[2] * rotationAxis[1] * (1 - cos(angle)) + rotationAxis[0] * sin(angle), cos(angle) + rotationAxis[2] * rotationAxis[2] * (1 - cos(angle));
    return rotationMatrix;
}

int main(int argc, char* argv[]) {

    gsInfo << "This is Cook's membrane benchmark with nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    real_t youngsModulus = 240.565e6;
    real_t poissonsRatio = 0.4;
    index_t numUniRef = 0;
    index_t numDegElev = 0;
    index_t numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
    cmd.addReal("p", "poisson", "Poisson's ratio used in the material law", poissonsRatio);
    cmd.addInt("r", "refine", "Number of uniform refinement application", numUniRef);
    cmd.addInt("d", "degelev", "Number of degree elevation application", numDegElev);
    cmd.addInt("s", "point", "Number of points to plot to Paraview", numPlotPoints);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    gsKnotVector<> KV1(0, 1, 0, 5);
    gsKnotVector<> KV2(0, 1, 0, 5);
    gsMatrix<real_t> C(25, 3);
    gsMatrix<real_t> CTemp(25, 3);
    gsMatrix<real_t> W(25, 1);
	
    W <<
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

    CTemp <<
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
    C = CTemp;
    gsTensorNurbs<2, real_t> patch1(KV1, KV2, C, W);
	
    gsVector<real_t> rotAxis(3);
    rotAxis << sqrt(3.) / 3., sqrt(3.) / 3., sqrt(3.) / 3.;
	
    gsMatrix<real_t> rotMat = generateRotationMatrix(M_PI * 2. / 3., rotAxis);
    C = CTemp * rotMat.transpose();
    gsTensorNurbs<2, real_t> patch2(KV1, KV2, C, W);
    rotMat = generateRotationMatrix(-M_PI * 2. / 3., rotAxis);
    C = CTemp * rotMat.transpose();
    gsTensorNurbs<2, real_t> patch3(KV1, KV2, C, W);
    // scanning geometry
    gsMultiPatch<> geometry;
    geometry.addPatch(patch1);
    geometry.addPatch(patch2);
    geometry.addPatch(patch3);
	
    geometry.computeTopology();
    // creating bases
    gsMultiBasis<> basisDisplacement(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basisDisplacement.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basisDisplacement.uniformRefine();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // neumann BC
    gsConstantFunction<> f(0., 625e4, 2);

    // boundary conditions
    gsBoundaryConditions<> bcInfo;
    for (index_t d = 0; d < 2; ++d)
        bcInfo.addCondition(0, boundary::west, condition_type::dirichlet, nullptr, d);
    bcInfo.addCondition(0, boundary::east, condition_type::neumann, &f);

    // source function, rhs
    gsConstantFunction<> g(0., 0., 2);

    //=============================================//
                  // Solving //
    //=============================================//

    // creating assembler
    gsElasticityAssemblerElasticSurface<real_t> assembler(geometry, basisDisplacement, bcInfo, g);
    assembler.options().setReal("YoungsModulus", youngsModulus);
    assembler.options().setReal("PoissonsRatio", poissonsRatio);
    assembler.options().setInt("MaterialLaw", material_law::neo_hooke_ln);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsIterative<real_t> solver(assembler);
    solver.options().setInt("Verbosity", solver_verbosity::all);
    solver.options().setInt("Solver", linear_solver::LDLT);

    //=============================================//
                  // Output //
    //=============================================//

    // solution to the nonlinear problem as an isogeometric displacement field
    gsMultiPatch<> displacement;
    assembler.constructSolution(solver.solution(), solver.allFixedDofs(), displacement);
    gsPiecewiseFunction<> geoError;
    assembler.constructGeoCalc(geoError, 1.);

    if (numPlotPoints > 0)
    {
        // constructing an IGA field (geometry + solution)
        gsField<> displacementField(assembler.patches(), displacement);
        gsField<> geoErrorField(assembler.patches(), geoError, false);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string, const gsField<>*> fields;
        fields["Displacement"] = &displacementField;
        fields["GeoError"] = &geoErrorField;
        gsWriteParaviewMultiPhysics(fields, "3patchesSphereSurface", numPlotPoints);
        gsInfo << "Open \"3patchesSphereSurface.pvd\" in Paraview for visualization.\n";
    }

    return 0;
}