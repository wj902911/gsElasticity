/** @file gsElasticityFunctions.hpp

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

#include <gsElasticity/gsElasticityFunctionsExtension.h>
#include <gsElasticity/gsWriteParaviewMultiPhysicsExtension.h>
#include <gsCore/gsFuncData.h>
#include <gsAssembler/gsAssembler.h>

namespace gismo
{

template <class T>
void gsCauchyStressFunctionExtension<T>::nonLinearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
#define REFERENCE 0
    result.setZero(m_dim,outputCols(u.cols()));
    // evaluating the fields
    gsMapData<T> mdGeo(NEED_GRAD_TRANSFORM);
    mdGeo.points = u;
    m_geometry->patch(m_patch).computeMap(mdGeo);
    gsMapData<T> mdDisp(NEED_DERIV);
    mdDisp.points = u;
    m_displacement->patch(m_patch).computeMap(mdDisp);
    // define temporary matrices here for efficieny
    gsMatrix<T> I = gsMatrix<T>::Identity(m_dim,m_dim);
#if REFERENCE
    gsMatrix<T> S, P, F, C, E;
#else
    gsMatrix<T> S,sigma,F,C,E;
#endif
    // material parameters
    T YM = m_options.getReal("YoungsModulus");
    T PR = m_options.getReal("PoissonsRatio");
    T lambda = YM * PR / ( ( 1. + PR ) * ( 1. - 2. * PR ) );
    T mu     = YM / ( 2. * ( 1. + PR ) );

    for (index_t q = 0; q < u.cols(); ++q)
    {
        // deformation gradient F = I + gradU*gradGeo^-1
        //if (mdGeo.jacobian(q).determinant() <= 0)
        if (mdGeo.jacobian(q).determinant() < 0)
            gsInfo << "Invalid domain parametrization: J = " << mdGeo.jacobian(q).determinant() <<
                      " at point (" << u.col(q).transpose() << ") of patch " << m_patch << std::endl;
        if (abs(mdGeo.jacobian(q).determinant()) > 1e-20)
            F = I + mdDisp.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
        else
            F = I;
        T J = F.determinant();
        if (J <= 0)
            gsInfo << "Invalid displacement field: J = " << J <<
                      " at point (" << u.col(q).transpose() << ") of patch " << m_patch << std::endl;
        // Second Piola-Kirchhoff stress tensor
        if (material_law::law(m_options.getInt("MaterialLaw")) == material_law::saint_venant_kirchhoff)
        {
            // Green-Lagrange strain tensor E = 0.5(F^T*F-I)
            E = (F.transpose() * F - I)/2;
            S = lambda*E.trace()*I + 2*mu*E;
        }
        if (material_law::law(m_options.getInt("MaterialLaw")) == material_law::neo_hooke_ln)
        {
            // Right Cauchy Green strain, C = F'*F
            C = F.transpose() * F;
            S = (lambda*log(J)-mu)*(C.cramerInverse()) + mu*I;
        }
        if (material_law::law(m_options.getInt("MaterialLaw")) == material_law::neo_hooke_quad)
        {
            // Right Cauchy Green strain, C = F'*F
            C = F.transpose() * F;
            S = (lambda*(J*J-1)/2-mu)*(C.cramerInverse()) + mu*I;
        }
        // transformation to P
#if REFERENCE
        P = F * S;
#else
        sigma = F*S*F.transpose()/J;
#endif
        gsVector<T> unormal;
        outerNormal(mdGeo, q, m_bs, unormal);
#if REFERENCE
        saveStress(P, result, q, unormal);
#else
        saveStress(sigma, result, q, unormal);
#endif
        //gsInfo << u.col(q) << "\n\n";
    }
}

template <class T>
void gsCauchyStressFunctionExtension<T>::saveStress(const gsMatrix<T> & S, gsMatrix<T> & result, index_t q, gsVector<T> unormal) const
{
    result.col(q) = S * unormal.normalized();
#if 0
    gsInfo << S << "\n\n";
    gsInfo << unormal.normalized() << "\n\n";
	
    switch (result.rows())
    {
    case 2:
        gsInfo << sqrt(result(0, q) * result(0, q) + result(1, q) * result(1, q)) << "\n\n";
		break;
	case 3:
		gsInfo << sqrt(result(0, q) * result(0, q) + result(1, q) * result(1, q) + result(2, q) * result(2, q)) << "\n\n";
    }
#endif
    
}

} // namespace gismo ends
