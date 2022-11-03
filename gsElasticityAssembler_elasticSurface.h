#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

namespace gismo
{
	
template <class T>
class gsElasticityAssembler_elasticSurface : public gsElasticityAssembler<T>
{
		
public:
	gsElasticityAssembler_elasticSurface(const gsMultiPatch<T>& patches,
			                             const gsMultiBasis<T>& basis,
			                             const gsBoundaryConditions<T>& bconditions,
			                             const gsFunction<T>& body_force);
		
	
		
protected:
	virtual void assemble(const gsMultiPatch<T>& displacement);
		
};

}


