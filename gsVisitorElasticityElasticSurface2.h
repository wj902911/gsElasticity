#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>
#include <fstream>

namespace gismo
{
	
template <class T>
class gsVisitorElasticityElasticSurface2
{
public:
    gsVisitorElasticityElasticSurface2(const gsPde<T>& pde_, 
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
        T YM = options.getReal("SurfaceYoungsModulus");
        T PR = options.getReal("SurfacePoissonsRatio");
        lambda = YM * PR / ((1. + PR) * (1. - 2. * PR));
        mu = YM / (2. * (1. + PR));
        gamma = options.getReal("SurfaceTension");
        I = gsMatrix<T>::Identity(dim, dim);
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
        short_t fixDir = patchSide.direction();
        // loop over quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {

            G_sub = md.jacobian(q);
            G_sub.removeCol(fixDir);
            g_sub = mdDisplacement.jacobian(q);
            g_sub.removeCol(fixDir);
            g_sub += G_sub;
            //calculate G_sub_alphaBeta and g_sub_alphaBeta.
            G_sub_ab = G_sub.transpose() * G_sub;
            g_sub_ab = g_sub.transpose() * g_sub;
            //calculate G_sup_alphaBeta and g_sup_alphaBeta.
            G_sup_ab = G_sub_ab.cramerInverse();
            g_sup_ab = g_sub_ab.cramerInverse();
            //calculate G_sup & g_sup, G_sup_1 & G_sup_2, g_sup_1 & g_sup_2 save in different columns in a single matrix.
            G_sup = G_sub * G_sup_ab;
            g_sup = g_sub * g_sup_ab;
            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            const index_t numGrads = basisValuesDisp[1].rows() / md.dim.first;
            const gsAsConstMatrix<T> grads_kTemp(basisValuesDisp[1].col(q).data(), md.dim.first, numGrads);
            gsMatrix<T> grads_k = grads_kTemp.transpose();
            grads_k.removeCol(fixDir);
            physGrad.noalias() = G_sup * grads_k.transpose();
            // deformation gradient F, f
            F = g_sub * G_sup.transpose();
            f = G_sub * g_sup.transpose();
            //calculate Jacobian J.
			switch(dim)
            {
            case 2:
                J = g_sub.norm() / G_sub.norm();
				break;
            case 3:
                const gsMatrix<T, 3, 1> g1 = g_sub.col(0);
                const gsMatrix<T, 3, 1> g2 = g_sub.col(1);
                const gsMatrix<T, 3, 1> G1 = G_sub.col(0);
                const gsMatrix<T, 3, 1> G2 = G_sub.col(1);
				J = g1.cross(g2).norm() / G1.cross(G2).norm();
                break;
            }
            // Right Cauchy Green strain, C = F'*F
            //RCG = F.transpose() * F;
            // Inverse of Right Cauchy Green strain, C = f*f'
            RCGinv = f * f.transpose();
            outerNormal(md, q, patchSide, unormal);
            const T weight = quWeights[q] * unormal.norm();
			//calculate Stress
            S = (lambda * log(J) - mu + gamma * J) * RCGinv + mu * I;
            // elasticity tensor
            matrixTraceTensor<T>(C, RCGinv, RCGinv);
            C *= lambda + gamma * J;
            symmetricIdentityTensor<T>(Ctemp, RCGinv);
            C += (mu - lambda * log(J) - gamma * J) * Ctemp;
            // loop over active basis functions (u_i)
            for (index_t i = 0; i < N_D; i++)
            {
                // Material tangent K_tg_mat = B_i^T * C * B_j;
                setB<T>(B_i, F, physGrad.col(i));
                materialTangentTemp = B_i.transpose() * C;
                // Geometric tangent K_tg_geo = gradB_i^T * S * gradB_j;
                geometricTangentTemp = S * physGrad.col(i);
                // loop over active basis functions (v_j)
                for (index_t j = 0; j < N_D; j++)
                {
                    setB<T>(B_j, F, physGrad.col(j));
                    materialTangent = materialTangentTemp * B_j;
                    T geometricTangent = geometricTangentTemp.transpose() * physGrad.col(j);
                    // K_tg = K_tg_mat + I*K_tg_geo;
                    for (short_t d = 0; d < dim; ++d)
                        materialTangent(d, d) += geometricTangent;
                    for (short_t di = 0; di < dim; ++di)
                        for (short_t dj = 0; dj < dim; ++dj)
                            localMat(di * N_D + i, dj * N_D + j) += weight * materialTangent(di, dj);
                }
                // Second Piola-Kirchhoff stress tensor as vector
                voigtStress<T>(Svec, S);
                // rhs = -r = force - B*Svec,
                localResidual = B_i.transpose() * Svec;
                for (short_t d = 0; d < dim; d++)
                    localRhs(d * N_D + i) -= weight * localResidual(d);
            }
        }
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
    }

protected:
    // problem info
    short_t dim;
    index_t patch; // current patch
    const gsBasePde<T>* pde_ptr;
    T lambda, mu, gamma, J;
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
    gsMatrix<T> I, C, Ctemp, physGrad, F, f, RCGinv, B_i, materialTangentTemp, B_j, materialTangent, S, G_sub, g_sub, G_sup, g_sup, G_sub_ab, g_sub_ab, G_sup_ab, g_sup_ab, S_sup_ab;
    gsVector<T> unormal, geometricTangentTemp, Svec, localResidual;
	
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;
};

}