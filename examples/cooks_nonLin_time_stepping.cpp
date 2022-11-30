/// This is the "Cook's membrane" benchmark solved using the nonlinear elasticity solver.
/// The problem description and reference solutions can be found in the Ph.D. thesis of O.Weeger
/// "Isogeometric Finite Element Analysis of Nonlinear Structural Vibrations", 2015.
///
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsIncrementIterative.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

using namespace gismo;

int main(int argc, char* argv[]){

    gsInfo << "This is Cook's membrane benchmark with nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = ELAST_DATA_DIR"/cooks.xml";
    real_t youngsModulus = 240.565e6;
    real_t poissonsRatio = 0.4;
    index_t numUniRef = 4;
    index_t numDegElev = 1;
    index_t numPlotPoints = 10000;

    // minimalistic user interface for terminal
    gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
    cmd.addReal("p","poisson","Poisson's ratio used in the material law",poissonsRatio);
    cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    cmd.addInt("s","point","Number of points to plot to Paraview",numPlotPoints);
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

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    real_t endTime = 1;
    int numSteps = 10;
    real_t Dt = endTime / numSteps ;
    gsMatrix<> totalsolVector;
    for ( int i = 0; i<=numSteps; ++i){
    // neumann BC
        gsConstantFunction<> fy(i*Dt*0.001,1);

        // boundary conditions
        gsBoundaryConditions<> bcInfo;
        for (index_t d = 0; d < 2; ++d){
            bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,d);
        }
        bcInfo.addCondition(0,boundary::east,condition_type::dirichlet,&fy);


        // source function, rhs
        gsConstantFunction<> g(0.,0.,0);

        //=============================================//
                      // Solving //
        //=============================================//

        // creating assembler
        gsElasticityAssembler<real_t> assembler(geometry,basisDisplacement,bcInfo,g);
        assembler.options().setReal("YoungsModulus",youngsModulus);
        assembler.options().setReal("PoissonsRatio",poissonsRatio);
        assembler.options().setInt("MaterialLaw",material_law::neo_hooke_ln);
        gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";
        totalsolVector.setZero(assembler.numDofs(),1);
        // setting Newton's method
        gsIncrementIterative<real_t> solver(assembler,totalsolVector);
        solver.options().setInt("Verbosity",solver_verbosity::all);
        solver.options().setInt("Solver",linear_solver::LDLT);

        gsInfo << "Solving...\n";
        gsStopwatch clock;
        clock.restart();
        solver.solve(totalsolVector);
        gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

        gsInfo<<"Norms is "<<totalsolVector.norm()<<"\n";

        //=============================================//
                      // Output //
        //=============================================//

        // solution to the nonlinear problem as an isogeometric displacement field
        gsMultiPatch<> displacement;
        assembler.constructSolution(totalsolVector,solver.allFixedDofs(),displacement);
        gsPiecewiseFunction<> stresses;
        assembler.constructCauchyStresses(displacement,stresses,stress_components::von_mises);

        if (numPlotPoints > 0)
        {
            // constructing an IGA field (geometry + solution)
            gsField<> displacementField(assembler.patches(),displacement);
            gsField<> stressField(assembler.patches(),stresses,true);
            // creating a container to plot all fields to one Paraview file
            std::map<std::string,const gsField<> *> fields;
            fields["Displacement"] = &displacementField;
            fields["von Mises"] = &stressField;
            gsWriteParaviewMultiPhysics(fields,"cooks"+std::to_string(i),numPlotPoints);
            gsInfo << "Open \"cooks.pvd\" in Paraview for visualization.\n";
        }
}

    return 0;
}
