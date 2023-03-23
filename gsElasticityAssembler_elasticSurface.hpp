#pragma once

#include <gsElasticity/gsElasticityAssembler_elasticSurface.h>

#include <gsElasticity/gsVisitorNonLinearElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>
#include <gsElasticity/gsVisitorElasticityElasticSurface.h>
#include <gsElasticity/gsVisitorElasticityElasticSurface2.h>


namespace gismo
{
	
template <class T>
gsElasticityAssemblerElasticSurface<T>::gsElasticityAssemblerElasticSurface(const gsMultiPatch<T>& patches,
		                                                                      const gsMultiBasis<T>& basis,
		                                                                      const gsBoundaryConditions<T>& bconditions,
		                                                                      const gsFunction<T>& body_force)
:gsElasticityAssembler<T>(patches, basis, bconditions, body_force)
{
    m_options.addReal("SurfaceTension", "Surface constant tension", 1.);
    m_options.addReal("SurfaceYoungsModulus", "Surface constant tension", 1.);
    m_options.addReal("SurfacePoissonsRatio", "Surface constant tension", .4);
}

template <class T>
void gsElasticityAssemblerElasticSurface<T>::constructCauchyStressesExtension(const gsMultiPatch<T>& displacement,
    gsPiecewiseFunction<T>& result,
    boundary::side s,
    stress_components::components comp) const
{
    if (comp == stress_components::all_2D_vector || comp == stress_components::all_2D_matrix)
        GISMO_ENSURE(m_dim == 2, "Invalid stress components for a 2D problem");
    if (comp == stress_components::normal_3D_vector || comp == stress_components::shear_3D_vector ||
        comp == stress_components::all_3D_matrix)
        GISMO_ENSURE(m_dim == 3, "Invalid stress type for a 3D problem");
    GISMO_ENSURE(m_options.getInt("MaterialLaw") == material_law::hooke ||
        m_options.getInt("MaterialLaw") == material_law::neo_hooke_ln ||
        m_options.getInt("MaterialLaw") == material_law::saint_venant_kirchhoff ||
        m_options.getInt("MaterialLaw") == material_law::neo_hooke_quad,
        "Pressure field not provided! Can't compute stresses with the chosen material law.");
    result.clear();

    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
        result.addPiecePointer(new gsCauchyStressFunctionExtension<T>(p, s, comp, m_options,
            &(m_pde_ptr->domain()), &displacement));
}

template<class T>
void gsElasticityAssemblerElasticSurface<T>::assemble(const gsMultiPatch<T>& displacement)
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
        gsVisitorElasticityElasticSurface2<T> surfaceVisitor(*m_pde_ptr, displacement, BCs[i]);
        Base::template push<gsVisitorElasticityElasticSurface2<T> >(surfaceVisitor, BCs[i]);
    }

    m_system.matrix().makeCompressed();

    //gsInfo << gsMatrix<T>(m_system.matrix()) << std::endl;
    //gsInfo << std::endl;

}

}
