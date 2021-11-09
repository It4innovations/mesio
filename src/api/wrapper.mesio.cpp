
#include "mesio.h"

#include "basis/containers/serializededata.h"
#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/utilities/sysutils.h"

#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.h"
#include "esinfo/config.h"
#include "esinfo/meshinfo.h"

#include "mesh/preprocessing/meshpreprocessing.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/bodystore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "wrappers/mpi/communication.h"

#include <cstdio>

struct MESIOData {
	MESIOData()
	{
		etypes.reserve(mesio::info::mesh->elements->distribution.process.size);
		for (size_t e = 0; e < mesio::info::mesh->elements->epointers->datatarray().size(); ++e) {
			etypes.push_back(static_cast<int>(mesio::info::mesh->elements->epointers->datatarray()[e]->code));
		}
		domains.reserve(mesio::info::mesh->elements->distribution.process.size);
		for (size_t i = 1; i < mesio::info::mesh->domains->elements.size(); ++i) {
			domains.resize(mesio::info::mesh->domains->elements[i], mesio::info::mesh->domains->offset + i - 1);
		}
		rtypes.resize(mesio::info::mesh->boundaryRegions.size());
		for (size_t i = 0; i < mesio::info::mesh->boundaryRegions.size(); ++i) {
			if (mesio::info::mesh->boundaryRegions[i]->dimension) {
				rtypes[i].reserve(mesio::info::mesh->boundaryRegions[i]->epointers->datatarray().size());
				for (size_t e = 0; e < mesio::info::mesh->boundaryRegions[i]->epointers->datatarray().size(); ++e) {
					rtypes[i].push_back(static_cast<int>(mesio::info::mesh->boundaryRegions[i]->epointers->datatarray()[e]->code));
				}
			}
		}
	}

	std::vector<esint> etypes, domains;
	std::vector<std::vector<esint> > rtypes;
};

using namespace mesio;

void MESIOInit(
	MPI_Comm		comm,
	int				verbosity)
{
	info::env::set();
	info::mpi::init(comm);
	MPITools::init();
	eslog::init(new Logger<ProgressTerminalLogger>);
	for (int i = 0; i < eslog::logger->size; ++i) {
		eslog::logger->args[i]->verbosity = verbosity;
	}
	Mesh::init();

	eslog::startln("MESIO: INITIALIZED", "MESIO");
}

void MESIOFinalize()
{
	Mesh::finish();
	MPITools::finish();
	eslog::endln("MESIO: FINISHED");
}

void MESIOLoad(
	MESIO*			mesio,
	MESIOFormat		format,
	const char*		path,
	MESIODecomposer	decomposer,
	int				domains)
{
	switch (format) {
	case MESIO_ANSYS: info::config::input.format = InputConfiguration::FORMAT::ANSYS_CDB; break;
	case MESIO_ENSIGHT: info::config::input.format = InputConfiguration::FORMAT::ENSIGHT; break;
	case MESIO_VTK_LEGACY: info::config::input.format = InputConfiguration::FORMAT::VTK_LEGACY; break;
	case MESIO_XDMF: info::config::input.format = InputConfiguration::FORMAT::XDMF; break;
	}
	info::config::input.path = path;
	switch (decomposer) {
	case MESIO_NONE: info::config::input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::NONE; break;
	case MESIO_METIS: info::config::input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::METIS; break;
	case MESIO_PARMETIS: info::config::input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::PARMETIS; break;
	case MESIO_PTSCOTCH: info::config::input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::PTSCOTCH; break;
	case MESIO_HILBERT_CURVE: info::config::input.decomposition.parallel_decomposer = DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE; break;
	}
	info::config::input.decomposition.domains = domains;
	info::mesh->preferedDomains = domains;

	Mesh::load();
	info::mesh->preprocess();

	*mesio = new MESIOData();
}

void MESIONodes(
	MESIO           mesio,
	MESIOInt*       nhalo,
	MESIOInt*       offset,
	MESIOInt*       size,
	MESIOInt*       totalSize,
	MESIOInt**      ids,
	MESIOInt**      position,
	MESIOReal**     coordinates)
{
	*nhalo = info::mesh->nodes->uniqInfo.nhalo;
	*offset = info::mesh->nodes->uniqInfo.offset;
	*size = info::mesh->nodes->uniqInfo.size;
	*totalSize = info::mesh->nodes->uniqInfo.totalSize;

	*ids = info::mesh->nodes->IDs->datatarray().data();
	*position = info::mesh->nodes->uniqInfo.position.data();
	*coordinates = static_cast<double*>(&info::mesh->nodes->coordinates->datatarray()[0].x);
}

void MESIONodesRanks(
	MESIO           mesio,
	MESIOInt**      rankDistribution,
	int**           rankData)
{
	*rankDistribution = info::mesh->nodes->ranks->boundarytarray().data();
	*rankData = info::mesh->nodes->ranks->datatarray().data();
}

void MESIONodesDomains(
	MESIO           mesio,
	MESIOInt**      domainDistribution,
	int**           domainData)
{
	*domainDistribution = info::mesh->nodes->domains->boundarytarray().data();
	*domainData = info::mesh->nodes->domains->datatarray().data();
}

void MESIONodeToElements(
	MESIO           mesio,
	MESIOInt**      elementDistribution,
	MESIOInt**      elementData
)
{

}

void MESIOElements(
	MESIO           mesio,
	MESIOInt*       offset,
	MESIOInt*       size,
	MESIOInt*       totalSize,
	MESIOInt**      type,
	MESIOInt**      enodesDistribution,
	MESIOInt**      enodesData)
{
	*offset = info::mesh->elements->distribution.process.offset;
	*size = info::mesh->elements->distribution.process.size;
	*totalSize = info::mesh->elements->distribution.process.totalSize;
	*type = mesio->etypes.data();
	*enodesDistribution = info::mesh->elements->nodes->boundarytarray().data();
	*enodesData = info::mesh->elements->nodes->datatarray().data();
}

void MESIOElementsDomains(
	MESIO           mesio,
	MESIOInt**      domains)
{
	*domains = mesio->domains.data();
}

void MESIOElementsMaterials(
	MESIO           mesio,
	int**           material)
{
	*material = info::mesh->elements->material->datatarray().data();
}

void MESIOElementsBodies(
	MESIO           mesio,
	int*            bodies,
	int**           body)
{
	*bodies = info::mesh->bodies->totalSize;
	*body = info::mesh->elements->body->datatarray().data();
}

void MESIOElementsFaceNeighbors(
	MESIO           mesio,
	MESIOInt**      neighborDistribution,
	MESIOInt**      neighborData)
{
	*neighborDistribution = info::mesh->elements->faceNeighbors->boundarytarray().data();
	*neighborData = info::mesh->elements->faceNeighbors->datatarray().data();
}

void MESIOElementsEdgeNeighbors(
	MESIO           mesio,
	MESIOInt**      neighborDistribution,
	MESIOInt**      neighborData
)
{
	*neighborDistribution = info::mesh->elements->edgeNeighbors->boundarytarray().data();
	*neighborData = info::mesh->elements->edgeNeighbors->datatarray().data();
}

void MESIOElementFaceList(
	MESIO               mesio,
	MESIOElementType    element,
	MESIOInt**          faceDistribution,
	MESIOInt**          faceData
)
{
	*faceDistribution = info::mesh->edata[element].faceList->boundarytarray().data();
	*faceData = info::mesh->edata[element].faceList->datatarray().data();
}

void MESIOElementEdgeList(
	MESIO               mesio,
	MESIOElementType    element,
	MESIOInt**          edgeDistribution,
	MESIOInt**          edgeData
)
{
	*edgeDistribution = info::mesh->edata[element].edgeList->boundarytarray().data();
	*edgeData = info::mesh->edata[element].edgeList->datatarray().data();
}

void MESIOElementsCounters(
	MESIO           mesio,
	MESIOInt        etype,
	MESIOInt*       offset,
	MESIOInt*       totalSize)
{
	*offset = info::mesh->elementsRegions.front()->distribution.code[etype].offset;
	*totalSize = info::mesh->elementsRegions.front()->distribution.code[etype].totalSize;
}

int MESIOElementsRegions(
	MESIO           mesio)
{
	return info::mesh->elementsRegions.size();
}

void MESIOElementsRegion(
	MESIO           mesio,
	MESIOInt        region,
	const char**    name,
	MESIOInt*       size,
	MESIOInt**      elements)
{
	*name = info::mesh->elementsRegions[region]->name.c_str();
	*size = info::mesh->elementsRegions[region]->elements->datatarray().size();
	*elements = info::mesh->elementsRegions[region]->elements->datatarray().data();
}

void MESIOElementsRegionNodes(
	MESIO           mesio,
	MESIOInt        region,
	MESIOInt*       nhalo,
	MESIOInt*       offset,
	MESIOInt*       size,
	MESIOInt*       totalSize,
	MESIOInt**      nodes,
	MESIOInt**      position)
{
	*nhalo = info::mesh->elementsRegions[region]->nodeInfo.nhalo;
	*offset = info::mesh->elementsRegions[region]->nodeInfo.offset;
	*size = info::mesh->elementsRegions[region]->nodeInfo.size;
	*totalSize = info::mesh->elementsRegions[region]->nodeInfo.totalSize;
	*nodes = info::mesh->elementsRegions[region]->nodes->datatarray().data();
	*position = info::mesh->elementsRegions[region]->nodeInfo.position.data();
}

void MESIOElementsRegionCounters(
	MESIO           mesio,
	MESIOInt        region,
	MESIOInt        etype,
	MESIOInt*       offset,
	MESIOInt*       totalSize)
{
	*offset = info::mesh->elementsRegions[region]->distribution.code[etype].offset;
	*totalSize = info::mesh->elementsRegions[region]->distribution.code[etype].totalSize;
}

int MESIOBoundaryRegions(
	MESIO           mesio)
{
	return info::mesh->boundaryRegions.size();
}

void MESIOBoundaryRegion(
	MESIO           mesio,
	MESIOInt        region,
	const char**    name,
	MESIOInt*       dimension,
	MESIOInt*       size,
	MESIOInt**      type,
	MESIOInt**      parent,
	MESIOInt**      elementDistribution,
	MESIOInt**      elementData)
{
	*name = info::mesh->boundaryRegions[region]->name.c_str();
	*dimension = info::mesh->boundaryRegions[region]->dimension;
	if (*dimension) {
		*size = info::mesh->boundaryRegions[region]->distribution.process.size;
		*type = mesio->rtypes[region].data();
		*parent = info::mesh->boundaryRegions[region]->emembership->datatarray().data();
		*elementDistribution = info::mesh->boundaryRegions[region]->elements->boundarytarray().data();
		*elementData = info::mesh->boundaryRegions[region]->elements->datatarray().data();
	} else {
		*size = 0;
		*type = *parent = *elementDistribution = *elementData = NULL;
	}
}

void MESIOBoundaryRegionNodes(
	MESIO           mesio,
	MESIOInt        region,
	MESIOInt*       nhalo,
	MESIOInt*       offset,
	MESIOInt*       size,
	MESIOInt*       totalSize,
	MESIOInt**      nodes,
	MESIOInt**      position)
{
	*nhalo = info::mesh->boundaryRegions[region]->nodeInfo.nhalo;
	*offset = info::mesh->boundaryRegions[region]->nodeInfo.offset;
	*size = info::mesh->boundaryRegions[region]->nodeInfo.size;
	*totalSize = info::mesh->boundaryRegions[region]->nodeInfo.totalSize;
	*nodes = info::mesh->boundaryRegions[region]->nodes->datatarray().data();
	*position = info::mesh->boundaryRegions[region]->nodeInfo.position.data();
}

void MESIOBoundaryRegionCounters(
	MESIO           mesio,
	MESIOInt        region,
	MESIOInt        etype,
	MESIOInt*       offset,
	MESIOInt*       totalSize)
{
	*offset = info::mesh->boundaryRegions[region]->distribution.code[etype].offset;
	*totalSize = info::mesh->boundaryRegions[region]->distribution.code[etype].totalSize;
}
