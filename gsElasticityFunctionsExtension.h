/** @file gsElasticityFunctions.h

    @brief Provides useful classes derived from gsFunction which can be used
    for visualization or coupling.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsCore/gsMultiPatch.h>
#include <gsElasticity/gsBaseUtils.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/** @brief Compute Cauchy stresses for a previously computed/defined displacement field.
 *         Can be pushed into gsPiecewiseFunction to construct gsField for visualization in Paraview.
*/
template <class T>
class gsCauchyStressFunctionExtension : public gsFunction<T>
{
public:

    gsCauchyStressFunctionExtension(index_t patch, boundary::side s,stress_components::components comp,
                           const gsOptionList & options,
                           const gsMultiPatch<T> * geometry,
                           const gsMultiPatch<T> * displacement)
        : m_geometry(geometry),
          m_displacement(displacement),
          m_patch(patch),
          m_dim(m_geometry->patch(m_patch).parDim()),
          m_options(options),
          m_type(comp),
          m_bs(s)
    {}



    virtual short_t domainDim() const
    {
        return m_geometry->patch(m_patch).parDim();
    }

    virtual short_t targetDim() const
    {
        switch(m_type)
        {
        case stress_components::von_mises: return 1;
        case stress_components::all_2D_vector: return 3;
        case stress_components::all_2D_matrix: return 2;
        case stress_components::normal_3D_vector: return 3;
        case stress_components::shear_3D_vector: return 3;
        case stress_components::all_3D_matrix: return 3;
        default: return 0;
        };
    }

    /** @brief Each column of the input matrix (u) corresponds to one evaluation point.
     *         Columns of the output matrix (result) correspond to a set of stress components for vector cases,
     *         or consist of one number in case of stress_type::von_mises,
     *         or form square matrices concatinated in the col-direction.
     */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        nonLinearElastic(u, result); 
        return;
    }

protected:

    /// size of the output matrix according to the m_type
    index_t outputCols(index_t inputCols) const
    {
        switch (m_type)
        {
        case stress_components::von_mises: return inputCols;
        case stress_components::all_2D_vector: return inputCols;
        case stress_components::all_2D_matrix: return 2*inputCols;
        case stress_components::normal_3D_vector: return inputCols;
        case stress_components::shear_3D_vector: return inputCols;
        case stress_components::all_3D_matrix: return 3*inputCols;
        default: return 0;
        }
    }
    /// save components of the stress tensor to the output matrix according to the m_type
    void saveStress(const gsMatrix<T> & S, gsMatrix<T> & result, index_t q, gsVector<T> unormal) const;

    /// computation routines for different material laws
    void nonLinearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const;

protected:
    const gsMultiPatch<T> * m_geometry;
    const gsMultiPatch<T> * m_displacement;
    index_t m_patch;
    short_t m_dim;
    const gsOptionList & m_options;
    stress_components::components m_type;
    boundary::side m_bs;
}; // class definition ends

template <class T>
class gsGeoCalcFunction : public gsFunction<T>
{
public:
	
    gsGeoCalcFunction(real_t radius)
        :m_radius(radius)
    {}

    virtual short_t domainDim() const { return 3; }

    virtual short_t targetDim() const { return 1; }
	
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result.setZero(targetDim(), u.cols());
        for (index_t q = 0; q < u.cols(); ++q)
        {
            result(0, q) = u.col(q).norm() - m_radius;
        }
    }
protected:
	real_t m_radius;
};


} // namespace ends


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityFunctions.hpp)
#endif
