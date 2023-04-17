/** @file gsWriteParaviewMultiPhysics.cpp

    @brief Provides implementation for gsWriteParaviewMultiPhysics.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy (TU Kaiserslautern)
    Inspired by gsWriteParaview.hpp by A. Mantzaflaris
*/

#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsWriteParaviewMultiPhysicsExtension.h>
#include <gsUtils/gsPointGrid.h>
#include <gsUtils/gsMesh/gsMesh.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsField.h>
#include <gsIO/gsWriteParaview.h>
#include <gsElasticity/gsGeoUtils.h>


#define PLOT_PRECISION 11


namespace gismo
{
template<class T>
void gsWriteParaviewMultiPhysicsTimeStepWithMesh(
    std::map<std::string, const gsField<T> *> fields,
    std::string const & fn,
    gsParaviewCollection& collection, 
    int time, 
    unsigned npts, 
    bool mesh)
{
    const unsigned numP = fields.begin()->second->patches().nPatches();
    std::string fileName = fn.substr(fn.find_last_of("/\\")+1); // file name without a path

    for ( size_t p = 0; p < numP; ++p)
    {
        gsWriteParaviewMultiPhysicsSinglePatch(fields,p,fn + util::to_string(time) + "_" + util::to_string(p),npts);
        collection.addTimestep(fileName + util::to_string(time),p,time,".vts");
    }
    if (mesh)
    {
        
		
        for (size_t p = numP; p < 2 * numP; ++p)
        {
            const gsBasis<>& dom = fields.begin()->second->isParametrized() ?
                fields.begin()->second->igaFunction(p - numP).basis() : fields.begin()->second->patch(p - numP).basis();

            gsWriteParaviewMultiPhysicsSingleMesh(fields, p - numP, fn + util::to_string(time) + "_" + util::to_string(p) + "_mesh", 4);
            collection.addTimestep(fileName + util::to_string(time), p, time, "_mesh.vtp");
        }
    }

}

	
template<class T>
void gsWriteParaviewMultiPhysicsSingleMesh(
    std::map<std::string, const gsField<T>*> fields,
    const unsigned patchNum,
    std::string const& fn,
    unsigned resolution)
{
    const gsBasis<>& basis = fields.begin()->second->isParametrized() ?
        fields.begin()->second->igaFunction(patchNum).basis() : fields.begin()->second->patch(patchNum).basis();
    const gsGeometry<T>& Geo = fields.begin()->second->patches().patch(patchNum);

    gsMesh<T> msh(basis, resolution);
    gsMatrix<> pts(Geo.domainDim(), msh.numVertices());
    const int n = pts.rows();
    for (int i = 0; i < msh.numVertices(); i++)
    {
        pts.col(i) = msh.vertex(i).topRows(n);
    }
    Geo.evaluateMesh(msh);
    gsMatrix<>eval_geo(Geo.domainDim(), msh.numVertices());
    for (int i=0;i< msh.numVertices();i++)
    {
        eval_geo.col(i) = msh.vertex(i).topRows(n);
    }
    std::string mfn(fn);
    mfn.append(".vtp");
    std::ofstream file(mfn.c_str());
    if (!file.is_open())
        gsWarn << "gsWriteParaview: Problem opening file \"" << fn << "\"" << std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision(PLOT_PRECISION);
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "<PolyData>\n";
    file << "<Piece NumberOfPoints=\"" << msh.numVertices() << "\" NumberOfVerts=\"0\" NumberOfLines=\""
        << msh.numEdges() << "\" NumberOfStrips=\"0\" NumberOfPolys=\"" << msh.numFaces() << "\">\n";

    std::map<std::string, gsMatrix<> > data;
    for (typename std::map<std::string, const gsField<T>*>::iterator it = fields.begin(); it != fields.end(); it++)
    {
        data[it->first] = it->second->isParametric() ?
            it->second->function(patchNum).eval(pts) : it->second->function(patchNum).eval(eval_geo);
        //gsInfo << eval_geo.transpose() << "\n\n";
        //gsInfo << data[it->first].transpose() << "\n\n";
        if (data[it->first].rows() == 2)
        {
            data[it->first].conservativeResize(3, eval_geo.cols());
            data[it->first].row(2).setZero();
        }
    }

    file << "<PointData>\n";
    for (typename std::map<std::string, gsMatrix<T> >::iterator it = data.begin(); it != data.end(); it++)
    {
        file << "<DataArray type=\"Float32\" Name=\"" << it->first << "\" format=\"ascii\" NumberOfComponents=\"" << (it->second.rows() == 1 ? 1 : 3) << "\">\n";
        if (it->second.rows() == 1)
            for (index_t j = 0; j < it->second.cols(); ++j)
                file << it->second.at(j) << " \n";
        else
        {
            for (index_t j = 0; j < it->second.cols(); ++j)
            {
                for (index_t i = 0; i != it->second.rows(); ++i)
                    file << it->second(i, j) << " ";
                for (index_t i = it->second.rows(); i < 3; ++i)
                    file << "0 ";
                file << "\n";
            }
        }
        file << "</DataArray>\n";
    }
    file << "</PointData>\n";
	
    file << "<Points>\n";
    file << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (index_t j = 0; j < eval_geo.cols(); ++j)
    {
        for (index_t i = 0; i != n; ++i)
            file << eval_geo(i, j) << " ";
        for (index_t i = n; i < 3; ++i)
            file << "0 ";
        file << "\n";
    }

    file << "\n";
    file << "</DataArray>\n";
    file << "</Points>\n";

    file << "<Lines>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (typename std::vector< gsEdge<T> >::const_iterator it = msh.edges().begin();
        it != msh.edges().end(); ++it)
    {
        file << it->source->getId() << " " << it->target->getId() << "\n";
    }
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int count = 0;
    for (typename std::vector< gsEdge<T> >::const_iterator it = msh.edges().begin();
        it != msh.edges().end(); ++it)
    {
        count += 2;
        file << count << " ";
    }
    file << "\n";
    file << "</DataArray>\n";
    file << "</Lines>\n";

    file << "<Polys>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (typename std::vector< gsFace<T>* >::const_iterator it = msh.faces().begin();
        it != msh.faces().end(); ++it)
    {
        for (typename std::vector< gsVertex<T>* >::const_iterator vit = (*it)->vertices.begin();
            vit != (*it)->vertices.end(); ++vit)
        {
            file << (*vit)->getId() << " ";
        }
        file << "\n";
    }
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    count = 0;
    for (typename std::vector< gsFace<T>* >::const_iterator it = msh.faces().begin();
        it != msh.faces().end(); ++it)
    {
        count += (*it)->vertices.size();
        file << count << " ";
    }
    file << "\n";
    file << "</DataArray>\n";
    file << "</Polys>\n";

    file << "</Piece>\n";
    file << "</PolyData>\n";
    file << "</VTKFile>\n";
    file.close();
}

template<class T>
void gsWriteHistoryOutputBoundaryResults(
    std::map<std::string, const gsField<T>*> fields,
    std::string fn,
    std::string fieldName,
    boundary::side s,
    unsigned npts)
{
    const unsigned numP = fields.begin()->second->patches().nPatches();
    T sumResult = 0;
    for (unsigned i = 0; i < numP; ++i)
    {
        sumResult += gsWriteHistoryOutputBoundaryResultsSinglePatch(fields, fieldName, s, i, npts);
    }
	std::ofstream of(fn.c_str(), std::ios::app);
    of << sumResult / numP << "\n";
    of.close();
}

template<class T>
T gsWriteHistoryOutputBoundaryResultsSinglePatch(
    std::map<std::string, const gsField<T>*> fields,
    std::string fieldName,
    boundary::side s,
    const unsigned patchNum, 
    unsigned npts)
{
    const gsGeometry<>& geometry = fields.begin()->second->patches().patch(patchNum);
    const short_t d = geometry.domainDim();
    gsMatrix<> ab = geometry.support();
    gsVector<> a = ab.col(0);
    gsVector<> b = ab.col(1);
    gsVector<unsigned> np = distributePoints<T>(geometry, npts);
    //gsInfo << ab << "\n\n";
    //gsInfo << np << "\n\n";
    if (d == 3)
    {
        switch (s)
        {
        case boundary::west:
            b(0) = 0;
            np(0) = 1;
            break;
        case boundary::east:
            a(0) = 1;
            np(0) = 1;
            break;
        case boundary::south:
            b(1) = 0;
            np(1) = 1;
            break;
        case boundary::north:
            a(1) = 1;
            np(1) = 1;
            break;
        case boundary::front:
            b(2) = 0;
            np(2) = 1;
            break;
        case boundary::back:
            a(2) = 1;
            np(2) = 1;
            break;
        }
    }
    else
    {
        switch (s)
        {
        case boundary::west:
            b(0) = 0;
            np(0) = 1;
            break;
        case boundary::east:
            a(0) = 1;
            np(0) = 1;
            break;
        case boundary::south:
            b(1) = 0;
            np(1) = 1;
            break;
        case boundary::north:
            a(1) = 1;
            np(1) = 1;
            break;
        }
    }
    //gsInfo << a << "\n\n";
    //gsInfo << b << "\n\n";
    //gsInfo << np << "\n\n";

    gsMatrix<> pts = gsPointGrid(a, b, np);
    gsMatrix<> eval_geo = geometry.eval(pts);
    //gsInfo << eval_geo << "\n\n";
    //gsInfo << eval_geo << "\n\n";
	gsMatrix<> outputData;
    for (typename std::map<std::string, const gsField<T>*>::iterator it = fields.begin(); it != fields.end(); it++)
    {
        if (it->first == fieldName)
        {
            outputData = it->second->isParametric() ?
                it->second->function(patchNum).eval(pts) : it->second->function(patchNum).eval(eval_geo);

            if (outputData.rows() == 2)
            {
                outputData.conservativeResize(3, eval_geo.cols());
                outputData.row(2).setZero();
            }
        }
    }
    gsMatrix<> totalDisp
        = sqrt(outputData.row(0).array() * outputData.row(0).array()
               + outputData.row(1).array() * outputData.row(1).array()
               + outputData.row(2).array() * outputData.row(2).array());
    return totalDisp.mean();
	
}

}
#undef PLOT_PRECISION
