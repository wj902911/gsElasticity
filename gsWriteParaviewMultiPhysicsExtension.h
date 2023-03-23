/** @file gsWriteParaviewMultiPhysics.h

    @brief Allows to write several fields defined on the same geometry
    in one file, making it easier to operate with them inside Paraview.
    Ideally should be a part of gismoIO module,

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy (TU Kaiserslautern)
    Inspired by gsWriteParaview.h by A. Mantzaflaris
*/
#pragma once

#include <gismo.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsIO/gsParaviewCollection.h>

#define NS 1000

namespace gismo
{
    /// \brief Write a file containing several fields defined on the same geometry to ONE paraview file
/// and adds it as a timestep to a Paraview collection
/// \param fields a map of field pointers
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each patch
/// \param mesh if true, the parameter mesh is plotted as well
    template<class T>
    void gsWriteParaviewMultiPhysicsTimeStepWithMesh(
        std::map<std::string, const gsField<T>*> fields,
        std::string const& fn,
        gsParaviewCollection& collection, 
        int time, 
        unsigned npts = NS, 
        bool mesh = false);
	
    /// \brief Write a file containing several fields defined on the same geometry to ONE paraview file
    ///
    /// \param fields a map of field pointers
    /// \param fn filename where paraview file is written
    /// \param npts number of points used for sampling each patch
    /// \param mesh if true, the parameter mesh is plotted as well
    template<class T>
    void gsWriteParaviewMultiPhysicsSingleMesh(
        std::map<std::string, const gsField<T>*> fields,
        const unsigned patchNum,
        std::string const& fn,
        unsigned resolution = 8);

    template<class T>
    void gsWriteHistoryOutputBoundaryResults(
        std::map<std::string, const gsField<T>*> fields,
        std::string fn,
        std::string fieldName,
        boundary::side s,
        unsigned npts);
	
    template<class T>
    T gsWriteHistoryOutputBoundaryResultsSinglePatch(
        std::map<std::string, const gsField<T>*> fields,
        std::string fieldName,
        boundary::side s,
        const unsigned patchNum,
        unsigned npts);
}
#undef NS
