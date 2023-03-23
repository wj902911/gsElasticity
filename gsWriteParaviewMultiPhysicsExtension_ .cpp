
#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsWriteParaviewMultiPhysicsExtension.h>
#include <gsElasticity/gsWriteParaviewMultiPhysicsExtension.hpp>

namespace gismo
{
TEMPLATE_INST
void gsWriteParaviewMultiPhysicsTimeStepWithMesh(
    std::map<std::string, const gsField<real_t> *> fields, 
    std::string const & fn,
    gsParaviewCollection & collection, 
    int time, 
    unsigned npts,
    bool mesh);

	
	
TEMPLATE_INST
void gsWriteParaviewMultiPhysicsSingleMesh(
    std::map<std::string, const gsField<real_t>*> fields,
    const unsigned patchNum,
    std::string const& fn,
    unsigned resolution);

TEMPLATE_INST
void gsWriteHistoryOutputBoundaryResults(
    std::map<std::string, const gsField<real_t>*> fields,
    std::string fn,
    std::string fieldName,
    boundary::side s,
    unsigned npts);

TEMPLATE_INST
real_t gsWriteHistoryOutputBoundaryResultsSinglePatch(
    std::map<std::string, const gsField<real_t>*> fields,
    std::string fieldName,
    boundary::side s,
    const unsigned patchNum,
    unsigned npts);
}
