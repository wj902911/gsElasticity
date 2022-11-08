#pragma once

#include <gsElasticity/gsElasticityAssembler_elasticSurface.h>

#include <gsElasticity/gsVisitorNonLinearElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>
#include <gsElasticity/gsVisitorElasticityElasticSurface.h>


namespace gismo
{
	
template <class T>
gsElasticityAssembler_elasticSurface<T>::gsElasticityAssembler_elasticSurface(const gsMultiPatch<T>& patches,
		                                                                      const gsMultiBasis<T>& basis,
		                                                                      const gsBoundaryConditions<T>& bconditions,
		                                                                      const gsFunction<T>& body_force)
:gsElasticityAssembler<T>(patches, basis, bconditions, body_force)
{
    m_options.addReal("SurfaceTension", "Surface constant tension", 1.);
}

template<class T>
void gsElasticityAssembler_elasticSurface<T>::assemble(const gsMultiPatch<T>& displacement)
{
    GISMO_ENSURE(m_options.getInt("MaterialLaw") == material_law::saint_venant_kirchhoff ||
        m_options.getInt("MaterialLaw") == material_law::neo_hooke_ln ||
        m_options.getInt("MaterialLaw") == material_law::neo_hooke_quad,
        "Material law not specified OR not supported!");
    m_system.matrix().setZero();
    reserve();
    m_system.rhs().setZero();

    // Compute volumetric integrals and write to the global linear system
    gsVisitorNonLinearElasticity<T> visitor(*m_pde_ptr, displacement);
    Base::template push<gsVisitorNonLinearElasticity<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    // change to reuse rhs from linear system
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());
    typename gsBoundaryConditions<T>::bcContainer BCs(m_pde_ptr->bc().robinSides());
    for (int i=0;i< BCs.size();i++)
    {
        gsVisitorElasticityElasticSurface<T> surfaceVisitor(*m_pde_ptr, displacement, BCs[i]);
        Base::template push<gsVisitorElasticityElasticSurface<T> >(surfaceVisitor, BCs[i]);
    }

    m_system.matrix().makeCompressed();
}

}
