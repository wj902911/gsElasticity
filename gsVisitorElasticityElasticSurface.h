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
        // resize containers for global indices
        globalIndices.resize(dim);
        blockNumbers.resize(dim);
    }

    inline void evaluate(const gsBasisRefs<T>& basisRefs,
                         const gsGeometry<T>& geo,
                         const gsMatrix<T>& quNodes)
    {
		
    }

    inline void assemble(gsDomainIterator<T>& element,
                         const gsVector<T>& quWeights)
    {
		
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >& eliminatedDofs,
                              gsSparseSystem<T>& system)
    {
		
    }

protected:
    // problem info
    short_t dim;
    index_t patch; // current patch
    const gsBasePde<T>* pde_ptr;
    T gamma;
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

    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;
};

}