#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{
	
template <class T>
class gsVisitorElasticityElasticSurface
{
public:
    gsVisitorElasticityElasticSurface(const gsPde<T>& pde_, 
                                      const gsMultiPatch<T>& displacement_,
                                      const boundary_condition<T>& s)
        : pde_ptr(static_cast<const gsBasePde<T>*>(&pde_)),
          displacement(displacement_),
		  patchSide(s.side()) {}
	
    void initialize(const gsBasisRefs<T>& basisRefs,
                    const index_t patchIndex,
                    const gsOptionList& options,
                    gsQuadRule<T>& rule)
    {
        // parametric dimension of the first displacement component
        dim = basisRefs.front().dim();
        // a quadrature rule is defined by the basis for the first displacement component.
        rule = gsQuadrature::get(basisRefs.front(), options, patchSide.direction());
        // saving necessary info
        patch = patchIndex;
        gamma = options.getReal("SurfaceTension");
        // resize containers for global indices
        globalIndices.resize(dim);
        blockNumbers.resize(dim);
    }

    inline void evaluate(const gsBasisRefs<T>& basisRefs,
                         const gsGeometry<T>& geo,
                         const gsMatrix<T>& quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_VALUE to get points in the physical domain for evaluation of the RHS
        // NEED_MEASURE to get the Jacobian determinant values for integration
        // NEED_GRAD_TRANSFORM to get the Jacobian matrix to transform gradient from the parametric to physical domain
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // find local indices of the displacement basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0), localIndicesDisp);
        N_D = localIndicesDisp.rows();
		
        // Evaluate displacement basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes, 1, basisValuesDisp);
        //gsInfo << basisValuesDisp[1] << std::endl;
        //gsInfo << std::endl;
        // Evaluate right-hand side at the image of the quadrature points
        //pde_ptr->rhs()->eval_into(md.values[0], forceValues);
        // store quadrature points of the element for displacement evaluation
        mdDisplacement.points = quNodes;
        // NEED_DERIV to compute deformation gradient
        mdDisplacement.flags = NEED_DERIV;
        // evaluate displacement gradient
        displacement.patch(patch).computeMap(mdDisplacement);
    }

    inline void assemble(gsDomainIterator<T>& element,
                         const gsVector<T>& quWeights)
    {
        // initialize matrices and rhs
        localMat.setZero(dim * N_D, dim * N_D);
        localRhs.setZero(dim * N_D, 1);
        G_sub_ab.setZero(dim - 1, dim - 1);
        g_sub_ab.setZero(dim - 1, dim - 1);
        G_sup.setZero(dim, dim - 1);
        g_sup.setZero(dim, dim - 1);
        A_sup_abcd.resize(dim - 1, dim - 1);
        for (int i = 0; i < dim - 1; i++)
        {
            for (int j = 0; j < dim - 1; j++)
            {
                A_sup_abcd(i, j).setZero(dim - 1, dim - 1);
            }
        }
        short_t fixDir = patchSide.direction();
        // loop over quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            //calculate G_sub & g_sub, G_sub_1 & G_sub_2, g_sub_1 & g_sub_2 save in different columns in a single matrix.
            //gsInfo << md.jacobian(q) << std::endl;
            //gsInfo << std::endl;
            G_sub = md.jacobian(q);
            G_sub.removeCol(fixDir);
            g_sub = mdDisplacement.jacobian(q);
            g_sub.removeCol(fixDir);

            /*
            if (fixDir < dim - 1)
            {
                G_sub.block(0, fixDir, dim, dim - 1 - fixDir) = G_sub.block(0, fixDir + 1, dim, dim - 1 - fixDir);
                g_sub.block(0, fixDir, dim, dim - 1 - fixDir) = g_sub.block(0, fixDir + 1, dim, dim - 1 - fixDir);
            }
            G_sub.conservativeResize(dim, dim - 1);
            g_sub.conservativeResize(dim, dim - 1);
            */
            
            
            g_sub += G_sub;
            //gsInfo << g_sub << std::endl;
            //gsInfo << std::endl;
            //calculate G_sub_alphaBeta and g_sub_alphaBeta.
            for (int i = 0; i < dim - 1; i++)
            {
                for (int j = 0; j < dim - 1; j++)
                {
                    G_sub_ab(i, j) = G_sub.col(i).transpose() * G_sub.col(j);
                    g_sub_ab(i, j) = g_sub.col(i).transpose() * g_sub.col(j);
                }
            }
            //calculate G_sup_alphaBeta and g_sup_alphaBeta.
            G_sup_ab = G_sub_ab.inverse();
            g_sup_ab = g_sub_ab.inverse();
            //calculate G_sup & g_sup, G_sup_1 & G_sup_2, g_sup_1 & g_sup_2 save in different columns in a single matrix.
            for (int i = 0; i < dim - 1; i++)
            {
                for (int j = 0; j < dim - 1; j++)
                {
                    G_sup.col(i) += G_sup_ab(i, j) * G_sub.col(i);
                    g_sup.col(i) += g_sup_ab(i, j) * g_sub.col(i);
                }
            }
			//calculate Jacobian J.
			switch(dim)
            {
            case 2:
                J = g_sub.norm() / G_sub.norm();
				break;
            case 3:
                //gsInfo << g_sub.col(0) << std::endl;
                //gsInfo << g_sub.col(1) << std::endl;
                const gsMatrix<T, 3, 1> g1 = g_sub.col(0);
                const gsMatrix<T, 3, 1> g2 = g_sub.col(1);
                const gsMatrix<T, 3, 1> G1 = G_sub.col(0);
                const gsMatrix<T, 3, 1> G2 = G_sub.col(1);
				J = g1.cross(g2).norm() / G1.cross(G2).norm();
                break;
            }
			//calculate Stress
            S_sup_ab = gamma * J * g_sup_ab;
            //gsInfo << S_sup_ab << std::endl;
            //gsInfo << std::endl;
            //calculate A
			for (int a = 0; a < dim - 1; a++)
			{
				for (int b = 0; b < dim - 1; b++)
				{
					for (int c = 0; c < dim - 1; c++)
					{
                        for (int d = 0; d < dim - 1; d++)
                        {
                            A_sup_abcd(a, b)(c, d) = gamma * J * (-g_sup_ab(a, c) * g_sup_ab(b, d) 
                                                                 - g_sup_ab(a, d) * g_sup_ab(b, c) 
                                                                 + g_sup_ab(a, b) * g_sup_ab(c, d));
                        }
					}
				}
			}
            outerNormal(md, q, patchSide, unormal);
            const T weight = quWeights[q] * unormal.norm();
            gsMatrix<T> der_ab = basisValuesDisp[1].reshapeCol(q, dim, basisValuesDisp[1].rows() / dim).transpose();
            der_ab.removeCol(fixDir);
            //gsInfo << der_ab << std::endl;
            //gsInfo << std::endl;
            //gsInfo << quWeights[q] << std::endl;
            //gsInfo << std::endl;
            //gsInfo << unormal.norm() << std::endl;
            //gsInfo << std::endl;
            //calculate local Mat
            for (index_t i = 0; i < N_D; i++)
            {
                if (der_ab(i, 0) == 0)
                {
					continue;
                }
                for (index_t j = 0; j < N_D; j++)
                {
                    gsMatrix<T> temp;
                    temp.setZero(dim, dim);
                    if (der_ab(j, 0) == 0)
                    {
                        continue;
                    }
                    for (int a = 0; a < dim - 1; a++)
                    {
                        for (int b = 0; b < dim - 1; b++)
                        {
                            for (int c = 0; c < dim - 1; c++)
                            {
                                for (int d = 0; d < dim - 1; d++)
                                {
                                    temp += A_sup_abcd(c, b)(a, d) * g_sub.col(c) * g_sub.col(d).transpose() * der_ab(j, a) * der_ab(i, b);
                                }
                            }
                        }
                    }
                    for (int a = 0; a < dim - 1; a++)
                    {
                        for (int b = 0; b < dim - 1; b++)
                        {
                            Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(dim, dim);
                            //gsInfo << S_sup_ab(a, b) * der_ab(j, a) * der_ab(i, b) * Id << std::endl;
                            temp += S_sup_ab(a, b) * der_ab(j, a) * der_ab(i, b) * Id;
                        }
                    }
                    //gsInfo << temp << std::endl;
                    //gsInfo << std::endl;
                    for (short_t di = 0; di < dim; ++di)
                        for (short_t dj = 0; dj < dim; ++dj)
                        {
                            localMat(di * N_D + i, dj * N_D + j) += weight * temp(di, dj);
                        }
                }
            }
            //calculate local RHS
            for (short_t d = 0; d < dim; ++d)
            {
                for (short_t a = 0; a < dim - 1; ++a)
                {
                    for (short_t b = 0; b < dim - 1; ++b)
                    {
                        localRhs.middleRows(d * N_D, N_D).noalias() -= weight * S_sup_ab(a, b) * g_sub(d, b) * der_ab.col(a);
                    }
                }
            }
            
            //gsInfo << localMat << std::endl;
            //gsInfo << std::endl;
            //gsInfo << localRhs << std::endl;
            //gsInfo << std::endl;
        }
        //gsInfo << localMat << std::endl;
        //gsInfo << std::endl;
        //gsInfo << localRhs << std::endl;
        //gsInfo << std::endl;
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >& eliminatedDofs,
                              gsSparseSystem<T>& system)
    {
        // computes global indices for displacement components
        for (short_t d = 0; d < dim; ++d)
        {
            system.mapColIndices(localIndicesDisp, patchIndex, globalIndices[d], d);
            blockNumbers.at(d) = d;
        }
        // push to global system
        system.pushToRhs(localRhs, globalIndices, blockNumbers);
        system.pushToMatrix(localMat, globalIndices, eliminatedDofs, blockNumbers, blockNumbers);

        //gsInfo << system.matrix() << std::endl;
        //gsInfo << std::endl;
        //gsInfo << system.rhs() << std::endl;
        //gsInfo << std::endl;
    }

protected:
    // problem info
    short_t dim;
    index_t patch; // current patch
    const gsBasePde<T>* pde_ptr;
    T gamma, J;
    //const gsFunction<T>* robinFunction_ptr;
    boxSide patchSide;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of the displacement basis functions active at the current element
    gsMatrix<index_t> localIndicesDisp;
    // number of displacement basis functions active at the current element
    index_t N_D;
    // values and derivatives of displacement basis functions at quadrature points at the current element
    // values are stored as a N_D x numQuadPoints matrix; not sure about derivatives, must be smth like N_D*dim x numQuadPoints
    std::vector<gsMatrix<T> > basisValuesDisp;

    // current displacement field
    const gsMultiPatch<T>& displacement;
    // evaluation data of the current displacement field
    gsMapData<T> mdDisplacement;
    // all temporary matrices defined here for efficiency
	gsMatrix<T> G_sub, g_sub, G_sup, g_sup, G_sub_ab, g_sub_ab, G_sup_ab, g_sup_ab, S_sup_ab;
    gsVector<T> unormal;
	gsMatrix<gsMatrix<T>> A_sup_abcd;
	
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;
};

}