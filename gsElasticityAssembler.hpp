/** @file gsElasticityAssembler.hpp

    @brief Provides nonlinear elasticity system matrices for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger
*/

#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>

// Element visitors
#include <gsElasticity/gsVisitorLinearElasticity.h>
#include <gsElasticity/gsVisitorNonLinElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>

// ---
#include <gsCore/gsBoxTopology.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>


namespace gismo
{

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(gsMultiPatch<T> const & patches,
												gsMultiBasis<T> const & bases,
												T E_modulus,
												T poissons_ratio,
												T density_rho,
												gsBoundaryConditions<T> const & bconditions,
												const gsFunction<T> & body_force
												// const gsMatrix<T> pointCoords,
												// const gsMatrix<T> pointForces
												)
:   Base(patches), 
	m_rho(density_rho),
    m_bConditions(bconditions),
    m_bodyForce(&body_force)
{
    // Copy the spline basis
    m_bases.push_back(bases);
    
    // Initialize material properties
    m_lambda = E_modulus * poissons_ratio / ( (1.+poissons_ratio)*(1.-2.*poissons_ratio)) ;
    m_mu     = E_modulus / (2.*(1.+poissons_ratio)) ;
    
	m_dim = body_force.targetDim();
	GISMO_ASSERT(m_dim == m_patches.parDim(), "Dimensions not matching!");
	GISMO_ASSERT(m_dim == m_bases.front().dim(), "Dimensions not matching!");
	GISMO_ASSERT(m_dim == 2 || m_dim == 3, "Only for 2D and 3D!");
	
    // Mark boundary dofs
    m_dofMappers.resize(m_dim);
	for (index_t i = 0; i < m_dim; i++)
		m_bases.front().getMapper(true, m_bConditions, m_dofMappers[i] ); 

    // Determine system size
    m_dofs = 0;
	for (index_t i = 0; i < m_dim; i++)
        m_dofs += m_dofMappers[i].freeSize();


    // Add shifts to global matrix indices to get the global dof
    // number for each unknown coordinate
	unsigned inShift = 0;
	unsigned bdShift = 0;

	for (index_t i = 0; i < m_dim; i++)
	{
        m_dofMappers[i].setShift( inShift );
		m_dofMappers[i].setBoundaryShift( bdShift );
		inShift += m_dofMappers[i].freeSize();
		bdShift += m_dofMappers[i].boundarySize();
	}
}

template<class T>
void gsElasticityAssembler<T>::assembleNeumann()
{
    std::cout << "Linear Elasticity: assemble Neumann BC." << std::endl;
	
	for ( typename gsBoundaryConditions<T>::const_iterator it = m_bConditions.neumannBegin();
          it != m_bConditions.neumannEnd(); 
		  ++it )
    {
        gsVisitorElasticityNeumann<T> neumann(*it->function(), it->side());

        // Note: it->unknown()
        this->apply(neumann, it->patch(), it->side() );
    }
}

template<class T>
void gsElasticityAssembler<T>::assemble()
{
    std::cout << "Linear Elasticity: assemble stiffness matrix." << std::endl;
	
	//computeDirichletDofs();
    index_t numDirichlet = 0;
    for (index_t i = 0; i < m_dim; ++i)
        numDirichlet += m_dofMappers[i].boundarySize();
    m_ddof.setZero(numDirichlet, 1);

    if (m_dofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    size_t nonZerosPerCol = m_dim;
    for (index_t i = 0; i < m_dim; ++i) // to do: improve
        nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;

    m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
    m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );
        
    // Resize the load vector
    m_rhs.setZero(m_dofs, 1 );

    // Assemble volume stiffness and load vector integrals
    gsVisitorLinearElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce);
    for (unsigned np=0; np < m_patches.nPatches(); ++np )
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Enforce Neumann boundary conditions
    assembleNeumann();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & deformed)
{
    std::cout << "Nonlinear elasticity: assemble residual and tangential stiffness matrix." << std::endl;
	
	if ( m_ddof.size() == 0 )
    {
        assemble();
        return;
    }

    // Initialize the matrix and rhs vector
    m_matrix.setZero();
    
    // Resize the load vector
    m_rhs.setZero(m_dofs, 1);

    gsVisitorNonLinElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce, deformed.patch(0));
	//gsVisitorNonLinElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce, deformed.patch(0));

    // Assemble volume stiffness and load vector integrals
    for (unsigned np=0; np < m_patches.nPatches(); ++np )
    {
        visitor.setDeformed( deformed.patch(np) );

        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_rhs
        this->apply(visitor, np);
    }

    // Enforce Neumann forces
    assembleNeumann();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();   
}

template<class T> // AM: is non-homogeneous not called yet
void gsElasticityAssembler<T>::computeDirichletDofsIntpl()
{
    const gsDofMapper &mapper = m_dofMappers.front();
    
    m_ddof.resize( mapper.boundarySize(), 1 ); //--mrhs
    
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd(); 
		  ++it )
    {
        const int unk = it->unknown();
        const int k   = it->patch();
        const gsBasis<T> & basis = (m_bases[unk])[k];

        // Get dofs on this boundary
        gsMatrix<unsigned> * boundary = basis.boundary(it->side()) ;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary->size(); ++i)
            {
                const int ii= mapper.bindex( (*boundary)(i) , k );
                m_ddof.row(ii).setZero();
            }
            delete boundary;
            continue;
        }

        // Get the side information
        int dir = it->side().direction( );
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( m_dim );

        for ( int i = 0; i < m_dim; ++i)
        {
            if ( i == dir )
            {
                gsVector<T> b(1); 
                b[0] = ( basis.component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {   
                rr.push_back( basis.component(i).anchors()->transpose() );
            }
        }

        // Compute dirichlet values
        gsMatrix<T> fpts = 
            it->function()->eval( m_patches[it->patch()].eval(  gsPointGrid( rr ) ) );

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = basis.boundaryBasis(it->side());
        gsGeometry<T> * geo = h->interpolate(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t k=0; k!= boundary->size(); ++k)
        {
            const int ii= mapper.bindex( (*boundary)(k) , it->patch() );
            m_ddof.row(ii) = dVals.row(k);
        }
        delete h;
        delete geo;
        delete boundary;
    }
}


template<class T>
void  gsElasticityAssembler<T>::setSolution(const gsMatrix<T>& solVector, 
                                            gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    std::vector<gsFunction<T> * > sols ;

    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();
		coeffs.setZero(); //OWEE

        for (index_t j = 0; j < m_dim; ++j)
        {
            const gsDofMapper & mapper = m_dofMappers[j];
            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
            }
        }
    }
}


template<class T>
void  gsElasticityAssembler<T>::updateSolution(const gsMatrix<T>& solVector, 
                                               gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    std::vector<gsFunction<T> * > sols ;

    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < m_dim; ++j)
        {
            const gsDofMapper & mapper = m_dofMappers[j];
            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
            }
        }
    }
}


template<class T>
void  gsElasticityAssembler<T>::constructSolution(const gsMatrix<T>& solVector, 
                                                  gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");
    
    // The final solution is the deformed shell, therefore we add the
    // solVector to the undeformed coefficients
    result = m_patches;

    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size(); // m_patches[p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < m_dim; ++j) // For all components x, y, z
        {
            const gsDofMapper & mapper = m_dofMappers[j];// grab mapper for this component

            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector ?
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    coeffs(i,j) = m_ddof(mapper.bindex(i, p), 0);
                }
            }
        }
    }
}


}// namespace gismo
