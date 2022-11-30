#pragma once
#include <gsElasticity/gsIterative.h>

namespace gismo
{

template <class T>
class gsBaseAssembler;

template <class T>
class gsIterative_multiStep : public gsIterative<T>
{
public:
	gsIterative_multiStep(gsBaseAssembler<T>& assembler_);
	
	void solve(index_t numSteps);

protected:
	index_t iStep;
};

}