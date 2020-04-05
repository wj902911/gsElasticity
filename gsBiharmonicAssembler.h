/** @file gsBiharmonicAssembler.h

    @brief Provides stiffness matrix for bi-harmonic equation in the mixed formulation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsBaseAssembler.h>

namespace gismo
{

template <class T>
class gsBiharmonicAssembler : public gsBaseAssembler<T>
{
public:
    typedef gsBaseAssembler<T> Base;

    /// @brief This assebmler uses mixed finite elements
    gsBiharmonicAssembler(const gsMultiPatch<T> & patches,
                          const gsMultiBasis<T> & basis,
                          const gsBoundaryConditions<T> & bconditions,
                          const gsFunction<T> & body_force);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// @brief Refresh routine to set dof-mappers
    virtual void refresh();

    /// @brief Assembles the matrix
    virtual void assemble();

    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDDoFs,
                          bool assembleMatrix = true) {assemble();}

    //--------------------- SOLUTION CONSTRUCTION ----------------------------------//

    /// @brief construct the solution of the equation
    virtual void constructSolution(const gsMatrix<T> & solVector,
                                   const std::vector<gsMatrix<T> > & fixedDoFs,
                                   gsMultiPatch<T> & solution) const;

    /// @brief construct the Laplacian of the solution
    virtual void constructSolutionAux(const gsMatrix<T> & solVector,
                                      const std::vector<gsMatrix<T> > & fixedDoFs,
                                      gsMultiPatch<T> & solutionAux) const;

    /// @brief construct both the solution and the Laplactian
    virtual void constructSolution(const gsMatrix<T> & solVector,
                                   const std::vector<gsMatrix<T> > & fixedDoFs,
                                   gsMultiPatch<T> & solutionMain, gsMultiPatch<T> & solutionAux) const;

protected:
    /// a custom reserve function to allocate memory for the sparse matrix
    virtual void reserve();

protected:
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_options;
    using Base::m_system;
    using Base::m_ddof;
};

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElMassAssembler.hpp)
#endif

