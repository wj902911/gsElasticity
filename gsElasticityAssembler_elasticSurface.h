#pragma once

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityFunctionsExtension.h>

namespace gismo
{
	
template <class T>
class gsElasticityAssemblerElasticSurface : public gsElasticityAssembler<T>
{
		
public:

	typedef gsBaseAssembler<T> Base;

	gsElasticityAssemblerElasticSurface(const gsMultiPatch<T>& patches,
			                             const gsMultiBasis<T>& basis,
			                             const gsBoundaryConditions<T>& bconditions,
			                             const gsFunction<T>& body_force);

	virtual void constructCauchyStressesExtension(const gsMultiPatch<T>& displacement,
		gsPiecewiseFunction<T>& result,
		boundary::side s,
		stress_components::components component = stress_components::von_mises) const;
		
protected:
	virtual void assemble(const gsMultiPatch<T>& displacement);

protected:
	using Base::m_options;

};

}


