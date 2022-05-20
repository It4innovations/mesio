
#include "mesh.h"

#include "esinfo/config.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/eslog.hpp"

#include "basis/containers/serializededata.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/packing.h"
#include "wrappers/mpi/communication.h"

#include "input/builders/input.h"
#include "input/parsers/ansyscdb/ansyscdb.h"
#include "input/parsers/openfoam/openfoam.h"
#include "input/parsers/abaqus/abaqus.h"
#include "input/parsers/xdmf/xdmf.h"
#include "input/parsers/ensight/ensight.h"
#include "input/parsers/vtklegacy/vtklegacy.h"
#include "input/parsers/netgen/netgen.h"
#include "input/parsers/neper/neper.h"

#include "preprocessing/meshpreprocessing.h"
#include "store/statisticsstore.h"
#include "store/elementstore.h"
#include "store/nodestore.h"
#include "store/domainstore.h"
#include "store/clusterstore.h"
#include "store/bodystore.h"
#include "store/elementsregionstore.h"
#include "store/boundaryregionstore.h"
#include "store/contactinterfacestore.h"
#include "store/surfacestore.h"
#include "store/contactstore.h"
#include "store/fetidatastore.h"

#include "output/output.h"
#include "output/visualization/debug.h"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

using namespace mesio;

Element Mesh::edata[(int)Element::CODE::SIZE];

bool Mesh::convertDatabase()
{
	return true;
}

void Mesh::init()
{
	edata[static_cast<int>(Element::CODE::POINT1   )]   .init<Element::CODE::POINT1   >();
	edata[static_cast<int>(Element::CODE::LINE2    )]   .init<Element::CODE::LINE2    >();
	edata[static_cast<int>(Element::CODE::TRIANGLE3)]   .init<Element::CODE::TRIANGLE3>();
	edata[static_cast<int>(Element::CODE::SQUARE4  )]   .init<Element::CODE::SQUARE4  >();
	edata[static_cast<int>(Element::CODE::TETRA4   )]   .init<Element::CODE::TETRA4   >();
	edata[static_cast<int>(Element::CODE::PYRAMID5 )]   .init<Element::CODE::PYRAMID5 >();
	edata[static_cast<int>(Element::CODE::PRISMA6  )]   .init<Element::CODE::PRISMA6  >();
	edata[static_cast<int>(Element::CODE::HEXA8    )]   .init<Element::CODE::HEXA8    >();
	edata[static_cast<int>(Element::CODE::LINE3    )]   .init<Element::CODE::LINE3    >();
	edata[static_cast<int>(Element::CODE::TRIANGLE6)]   .init<Element::CODE::TRIANGLE6>();
	edata[static_cast<int>(Element::CODE::SQUARE8  )]   .init<Element::CODE::SQUARE8  >();
	edata[static_cast<int>(Element::CODE::TETRA10  )]   .init<Element::CODE::TETRA10  >();
	edata[static_cast<int>(Element::CODE::PYRAMID13)]   .init<Element::CODE::PYRAMID13>();
	edata[static_cast<int>(Element::CODE::PRISMA15 )]   .init<Element::CODE::PRISMA15 >();
	edata[static_cast<int>(Element::CODE::HEXA20   )]   .init<Element::CODE::HEXA20   >();

	info::mesh = new Mesh();
}

void Mesh::load()
{
	MeshBuilder *data = NULL;
	switch (info::config::input.format) {
	case InputConfiguration::FORMAT::ANSYS_CDB:      data = new AnsysCDBLoader     (info::config::input); break;
	case InputConfiguration::FORMAT::OPENFOAM:       data = new OpenFOAMLoader     (info::config::input); break;
	case InputConfiguration::FORMAT::ABAQUS:         data = new AbaqusLoader       (info::config::input); break;
	case InputConfiguration::FORMAT::XDMF:           data = new XDMFLoader         (info::config::input); break;
	case InputConfiguration::FORMAT::ENSIGHT:        data = new EnsightLoader      (info::config::input); break;
	case InputConfiguration::FORMAT::VTK_LEGACY:     data = new VTKLegacyLoader    (info::config::input); break;
	case InputConfiguration::FORMAT::NETGET:         data = new NetgenNeutralLoader(info::config::input); break;
	case InputConfiguration::FORMAT::NEPER:          data = new NeperLoader        (info::config::input); break;
	}

	data->load();
	data->build();

	delete data;
}

void Mesh::finish()
{
	delete info::mesh;
}

Mesh::Mesh()
: elements(new ElementStore()), nodes(new NodeStore()),
  domains(new DomainStore()),
  clusters(new ClusterStore()),
  bodies(new BodyStore()),
  FETIData(new FETIDataStore()),
  surface(new SurfaceStore()), domainsSurface(new SurfaceStore()),
  contact(new ContactStore()),

  output(new Output()),
  _omitClusterization(false),
  _omitDecomposition(false)
{
	dimension = 3;
	preferedDomains = info::config::input.decomposition.domains;
	if (preferedDomains == 0) {
		preferedDomains = info::env::threads; // TODO: set better value;
	}
}

Mesh::~Mesh()
{
	// we need to delete output first in order to wait for unfinished output operations
	delete output;

	delete elements;
	delete nodes;
	delete domains;
	delete clusters;
	delete bodies;
	delete FETIData;
	delete surface;
	delete domainsSurface;
	delete contact;

	for (size_t i = 0; i < contactInterfaces.size(); ++i) {
		delete contactInterfaces[i];
	}
	for (size_t i = 0; i < boundaryRegions.size(); ++i) {
		delete boundaryRegions[i];
	}
	for (size_t i = 0; i < elementsRegions.size(); ++i) {
		delete elementsRegions[i];
	}
}

ElementsRegionStore* Mesh::eregion(const std::string &name)
{
	return elementsRegions[eregionIndex(name)];
}

BoundaryRegionStore* Mesh::bregion(const std::string &name)
{
	return boundaryRegions[bregionIndex(name)];
}

size_t Mesh::eregionIndex(const std::string &name)
{
	for (size_t r = 0; r < elementsRegions.size(); r++) {
		if (StringCompare::caseSensitiveEq(elementsRegions[r]->name, name)) {
			return r;
		}
	}
	eslog::error("Unknown region of elements with name '%s'.\n", name.c_str());
	return 0;
}

size_t Mesh::bregionIndex(const std::string &name)
{
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (StringCompare::caseSensitiveEq(boundaryRegions[r]->name, name)) {
			return r;
		}
	}
	eslog::error("Unknown boundary region '%s'\n.", name.c_str());
	return 0;
}

bool Mesh::onAllElements(const std::string &eregion) const
{
	return StringCompare::caseInsensitiveEq(eregion, "ALL_ELEMENTS");
}

bool Mesh::hasPhaseChange() const
{
//	for (size_t m = 0; m < materials.size(); m++) {
//		if (materials[m]->phase_change) {
//			return true;
//		}
//	}
	return false;
}

void Mesh::setMaterials()
{
//	materials.clear();
//	std::map<std::string, int> matindex;
//	for (auto mat = info::config::getPhysics()->materials.begin(); mat != info::config::getPhysics()->materials.end(); ++mat) {
//		mat->second.name = mat->first;
//		materials.push_back(&mat->second);
//		matindex[mat->first] = materials.size() - 1;
//	}
//
//	for (auto mat = info::config::getPhysics()->material_set.begin(); mat != info::config::getPhysics()->material_set.end(); ++mat) {
//		ElementsRegionStore *region = eregion(mat->first);
//		if (matindex.find(mat->second) == matindex.end()) {
//			eslog::globalerror("Unknown material '%s'.\n", mat->second.c_str());
//		}
//		int material = matindex.find(mat->second)->second;
//		for (auto e = region->elements->datatarray().cbegin(); e != region->elements->datatarray().cend(); ++e) {
//			elements->material->datatarray()[*e] = material;
//		}
//	}
}

void Mesh::analyze()
{
	// TODO: resolve the problem with non-continuity with more bodies
//	if (_withFETI) {
//		info::config::input.decomposition.force_continuity = true;
//	}

	_omitClusterization = false;
	_omitDecomposition = false;

	if (info::mpi::size == 1) {
		_omitClusterization = true;
	}
	if (info::config::input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::NONE) {
		_omitClusterization = true;
	}
	if (
			info::config::input.decomposition.separate_materials ||
			info::config::input.decomposition.separate_regions ||
			info::config::input.decomposition.separate_etypes) {
		_omitDecomposition = false;
	}
}

void Mesh::reclusterize()
{
	mesh::computeElementsFaceNeighbors(nodes, elements, neighbors);
//	if (_withEdgeDual) {
//		mesh::computeElementsEdgeNeighbors(nodes, elements, neighbors);
//	}

	if (!_omitClusterization) {
		if (info::config::input.decomposition.parallel_decomposer == DecompositionConfiguration::ParallelDecomposer::HILBERT_CURVE) {
			mesh::computeElementsCenters(nodes, elements);
		}
		std::vector<esint> partition;
		mesh::computeElementsClusterization(elements, nodes, partition);
		mesh::exchangeElements(elements, nodes, elementsRegions, boundaryRegions, neighbors, neighborsWithMe, partition);
		if (info::config::input.decomposition.force_continuity) {
			std::vector<int> component;
			esint csize = mesh::getStronglyConnectedComponents(elements, component);
			esint coffset = csize;
			esint clusters = Communication::exscan(coffset);
			if (clusters > info::mpi::size) {
				std::vector<esint> dualDist, dualData;
				mesh::computeComponentDual(elements, coffset, csize, component, neighbors, dualDist, dualData);
				mesh::computeContinuousClusterization(elements, nodes, dualDist, dualData, coffset, csize, component, neighborsWithMe, partition);
				mesh::exchangeElements(elements, nodes, elementsRegions, boundaryRegions, neighbors, neighborsWithMe, partition);
			}
		}
	}

	esint minsize = elements->IDs->datatarray().size();
	Communication::allReduce(&minsize, NULL, 1, MPITools::getType<esint>().mpitype, MPI_MIN);
	if (minsize == 0) {
		eslog::globalerror("MESIO quit computation: process without any elements detected.\n");
	}

	mesh::sortNodes(nodes, elements, boundaryRegions);

	// to be removed
//	NamedData *orientation = NULL, *poly = NULL;
//	if (info::config::input.insert_orientation) {
//		orientation = info::mesh->elements->appendData(3, NamedData::DataType::VECTOR, "ORIENTATION");
//		poly = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "POLY");
//	}
}

void Mesh::computePersistentParameters()
{
	setMaterials();
	mesh::fillRegionMask(elements, elementsRegions);
	mesh::processNamelessElements(elements, elementsRegions);
	mesh::computeElementDistribution(elements);
	mesh::computeRegionsElementDistribution(elements, elementsRegions);

	{ // compute boundary element from nodes
		ElementStore *halo = NULL;
		for (size_t r = 0; halo == NULL && r < boundaryRegions.size(); ++r) {
			if (boundaryRegions[r]->originalDimension < boundaryRegions[r]->dimension) {
				halo = mesh::exchangeHalo(elements, nodes, neighbors);
			}
		}
		for (size_t r = 0; r < boundaryRegions.size(); ++r) {
			if (boundaryRegions[r]->originalDimension < boundaryRegions[r]->dimension) {
				mesh::computeRegionsBoundaryElementsFromNodes(nodes, elements, halo, elementsRegions, boundaryRegions[r]);
			}
		}

		if (halo) {
			delete halo;
		}
	}

	mesh::computeBodies(elements, bodies, elementsRegions, neighbors);

//	if (info::config::input.contact_interfaces.size()) {
//		mesh::computeBodiesSurface(nodes, elements, elementsRegions, surface, neighbors);
//		mesh::computeWarpedNormals(surface);
//		mesh::exchangeContactHalo(surface, contact);
//		mesh::findCloseElements(contact);
//		mesh::computeContactInterface(surface, contact);
//		mesh::arrangeContactInterfaces(contact, bodies, elementsRegions, contactInterfaces);
//		eslog::checkpointln("MESH: CONTACT INTERFACE COMPUTED");
//	}

	mesh::computeRegionsBoundaryDistribution(nodes, boundaryRegions, contactInterfaces);

	mesh::computeRegionsElementNodes(nodes, elements, neighbors, elementsRegions);
	mesh::computeRegionsBoundaryNodes(neighbors, nodes, boundaryRegions, contactInterfaces);
	mesh::computeRegionsBoundaryParents(nodes, elements, boundaryRegions, contactInterfaces);

	if (dimension == 3 && info::config::output.format == OutputConfiguration::FORMAT::STL_SURFACE) {
		mesh::computeBodiesSurface(nodes, elements, elementsRegions, surface, neighbors);
		mesh::triangularizeSurface(surface);
		eslog::checkpointln("MESH: BODIES SURFACE COMPUTED");
	}

	eslog::checkpointln("MESH: BODIES COMPUTED");
	eslog::checkpointln("MESH: BOUNDARY REGIONS COMPOSED");
}

void Mesh::partitiate(int ndomains)
{
	std::vector<esint> permutation;
	std::vector<size_t> distribution;
	if (_omitDecomposition) {
		permutation.resize(elements->distribution.process.size);
		std::iota(permutation.begin(), permutation.end(), 0);
		esint psize = elements->distribution.process.size / ndomains;
		distribution.push_back(0);
		for (esint p = 0, offset = 0; p < ndomains; ++p, offset += psize) {
			distribution.push_back(offset + psize);
		}
		_omitDecomposition = false; // only the first run
	} else {
		mesh::computeElementsDecomposition(elements, ndomains, distribution, permutation);
	}
	mesh::permuteElements(elements, nodes, domains, elementsRegions, boundaryRegions, contactInterfaces, neighbors, distribution, permutation);
	mesh::computeDomainDual(nodes, elements, domains, neighbors, neighborsWithMe);
	mesh::computeClustersDistribution(domains, clusters);
	eslog::checkpointln("MESH: MESH DECOMPOSED");

	mesh::computeElementIntervals(domains, elements);
	mesh::computeRegionsElementIntervals(elements, elementsRegions);
	mesh::computeRegionsBoundaryIntervals(domains, boundaryRegions, contactInterfaces);
	eslog::checkpointln("MESH: ELEMENT REGIONS ARRANGED");
}

void Mesh::preprocess()
{
	analyze();

	eslog::startln("MESH: PREPROCESSING STARTED", "MESHING");

	reclusterize();
	computePersistentParameters();
	partitiate(preferedDomains);

	DebugOutput::mesh();
	eslog::endln("MESH: PREPROCESSING FINISHED");
}

void Mesh::duplicate()
{
	eslog::startln("MESH: CREATE DUPLICATED INSTANCES", "DUPLICATION");

	size_t packedSize = 0;

	if (info::mpi::irank == 0) {
		packedSize += utils::packedSize(dimension);
		packedSize += utils::packedSize(preferedDomains);

		packedSize += elements->packedFullSize();
		packedSize += nodes->packedFullSize();

		packedSize += utils::packedSize(elementsRegions.size());
		for (size_t i = 0; i < elementsRegions.size(); i++) {
			packedSize += elementsRegions[i]->packedFullSize();
		}
		packedSize += utils::packedSize(boundaryRegions.size());
		for (size_t i = 0; i < boundaryRegions.size(); i++) {
			packedSize += boundaryRegions[i]->packedFullSize();
		}
		packedSize += utils::packedSize(contactInterfaces.size());
		for (size_t i = 0; i < contactInterfaces.size(); i++) {
			packedSize += contactInterfaces[i]->packedFullSize();
		}

		packedSize += domains->packedFullSize();
		packedSize += clusters->packedFullSize();
		packedSize += bodies->packedFullSize();

		packedSize += FETIData->packedFullSize();

		packedSize += surface->packedFullSize();
		packedSize += domainsSurface->packedFullSize();
		packedSize += contact->packedFullSize();

		packedSize += utils::packedSize(neighbors);
		packedSize += utils::packedSize(neighborsWithMe);
		packedSize += utils::packedSize(_omitClusterization);
		packedSize += utils::packedSize(_omitDecomposition);
	}

	Communication::broadcast(&packedSize, sizeof(size_t), MPI_BYTE, 0, MPITools::instances);
	char *buffer = new char[packedSize];

	if (info::mpi::irank == 0) {
		char *p = buffer;
		utils::pack(dimension, p);
		utils::pack(preferedDomains, p);

		elements->packFull(p);
		nodes->packFull(p);

		utils::pack(elementsRegions.size(), p);
		for (size_t i = 0; i < elementsRegions.size(); i++) {
			elementsRegions[i]->packFull(p);
		}
		utils::pack(boundaryRegions.size(), p);
		for (size_t i = 0; i < boundaryRegions.size(); i++) {
			boundaryRegions[i]->packFull(p);
		}
		utils::pack(contactInterfaces.size(), p);
		for (size_t i = 0; i < contactInterfaces.size(); i++) {
			contactInterfaces[i]->packFull(p);
		}

		domains->packFull(p);
		clusters->packFull(p);
		bodies->packFull(p);

		FETIData->packFull(p);

		surface->packFull(p);
		domainsSurface->packFull(p);
		contact->packFull(p);

		utils::pack(neighbors, p);
		utils::pack(neighborsWithMe, p);
		utils::pack(_omitClusterization, p);
		utils::pack(_omitDecomposition, p);
	}

	eslog::checkpoint("MESH: MESH PACKED");
	eslog::param("size[MB]", packedSize);
	eslog::ln();

	Communication::broadcast(buffer, packedSize, MPI_CHAR, 0, MPITools::instances);

	eslog::checkpointln("MESH: PACKED DATA BROADCASTED");

	if (info::mpi::irank != 0) {
		for (size_t i = 0; i < elementsRegions.size(); i++) {
			delete elementsRegions[i];
		}
		elementsRegions.clear();
		for (size_t i = 0; i < boundaryRegions.size(); i++) {
			delete boundaryRegions[i];
		}
		boundaryRegions.clear();
		for (size_t i = 0; i < contactInterfaces.size(); i++) {
			delete contactInterfaces[i];
		}
		contactInterfaces.clear();

		const char *p = buffer;
		utils::unpack(dimension, p);
		utils::unpack(preferedDomains, p);

		elements->unpackFull(p);
		nodes->unpackFull(p);

		size_t size;
		utils::unpack(size, p);
		for (size_t i = 0; i < size; i++) {
			elementsRegions.push_back(new ElementsRegionStore(p));
		}
		utils::unpack(size, p);
		for (size_t i = 0; i < size; i++) {
			boundaryRegions.push_back(new BoundaryRegionStore(p));
		}
		utils::unpack(size, p);
		for (size_t i = 0; i < size; i++) {
			contactInterfaces.push_back(new ContactInterfaceStore(p));
		}

		domains->unpackFull(p);
		clusters->unpackFull(p);
		bodies->unpackFull(p);

		FETIData->unpackFull(p);

		surface->unpackFull(p);
		domainsSurface->unpackFull(p);
		contact->unpackFull(p);

		utils::unpack(neighbors, p);
		utils::unpack(neighborsWithMe, p);
		utils::unpack(_omitClusterization, p);
		utils::unpack(_omitDecomposition, p);
		setMaterials();
	}

	delete[] buffer;

	eslog::endln("MESH: DUPLICATION FINISHED");
}

void Mesh::toBuffer()
{
	for (size_t i = 0; i < elements->data.size(); ++i) {
		if (elements->data[i]->name.size()) {
			elements->data[i]->toBuffer();
		}
	}
	for (size_t i = 0; i < nodes->data.size(); ++i) {
		if (nodes->data[i]->name.size()) {
			nodes->data[i]->toBuffer();
		}
	}
}

void Mesh::printMeshStatistics()
{
	size_t namesize = 56;

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (namesize < boundaryRegions[r]->name.size() + 16) {
			namesize = boundaryRegions[r]->name.size() + 16;
		}
	}

	for (size_t r = 0; r < elementsRegions.size(); r++) {
		if (namesize < elementsRegions[r]->name.size() + 16) {
			namesize = elementsRegions[r]->name.size() + 16;
		}
	}

	auto ename = [] (int code) -> std::string {
		switch (static_cast<Element::CODE>(code)) {

		case Element::CODE::POINT1: return "POINT1";

		case Element::CODE::LINE2: return "LINE2";

		case Element::CODE::TRIANGLE3: return "TRIANGLE3";
		case Element::CODE::SQUARE4: return "SQUARE4";

		case Element::CODE::TETRA4: return "TETRA4";
		case Element::CODE::PYRAMID5: return "PYRAMID5";
		case Element::CODE::PRISMA6: return "PRISMA6";
		case Element::CODE::HEXA8: return "HEXA8";

		case Element::CODE::LINE3: return "LINE3";

		case Element::CODE::TRIANGLE6: return "TRIANGLE6";
		case Element::CODE::SQUARE8: return "SQUARE8";

		case Element::CODE::TETRA10: return "TETRA10";
		case Element::CODE::PYRAMID13: return "PYRAMID13";
		case Element::CODE::PRISMA15: return "PRISMA15";
		case Element::CODE::HEXA20: return "HEXA20";

		default:
			eslog::internalFailure("unknown element code.\n");
			return "";
		}
	};

	size_t nelements = 0, nelementsnodes = 0;
	size_t fregs = 0, nfaces = 0, nfacenodes = 0;
	size_t eregs = 0, nedges = 0, nedgenodes = 0;
	size_t nregs = 0, nnodes = 0;

	for (size_t r = 1; r < elementsRegions.size(); r++) {
		nelements += elementsRegions[r]->distribution.process.totalSize; nelementsnodes += elementsRegions[r]->nodeInfo.totalSize;
	}
	for (size_t r = 1; r < boundaryRegions.size(); r++) {
		switch (boundaryRegions[r]->originalDimension) {
		case 2: ++fregs; nfaces += boundaryRegions[r]->distribution.process.totalSize; nfacenodes += boundaryRegions[r]->nodeInfo.totalSize; break;
		case 1: ++eregs; nedges += boundaryRegions[r]->distribution.process.totalSize; nedgenodes += boundaryRegions[r]->nodeInfo.totalSize; break;
		case 0: ++nregs; nnodes += boundaryRegions[r]->nodeInfo.totalSize; break;
		default: break;
		}
	}

	esint ecountTotal = 0, scountTotal = 0;
	for (size_t b = 0; b < elementsRegions[0]->bodies.size(); ++b) {
		ecountTotal += elementsRegions[0]->bodyElements[b];
		scountTotal += elementsRegions[0]->bodyFaces[b];
	}

	switch (info::config::output.logger) {
	case OutputConfiguration::LOGGER::USER:
		eslog::info(" ====================================== MESH STATISTICS ====================================== \n");
		eslog::info(" ============================================================================================= \n");
		eslog::info("  %s%*s : %16s %16s %16s\n", "GEOMETRY", namesize - 26, " ", "MIN", "MAX", "LENGTH");
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		eslog::info("  %*s : %16f %16f %16f\n", namesize - 18, "X-COORDINATES", elementsRegions.front()->nodeInfo.min.x, elementsRegions.front()->nodeInfo.max.x, elementsRegions.front()->nodeInfo.max.x - elementsRegions.front()->nodeInfo.min.x);
		eslog::info("  %*s : %16f %16f %16f\n", namesize - 18, "Y-COORDINATES", elementsRegions.front()->nodeInfo.min.y, elementsRegions.front()->nodeInfo.max.y, elementsRegions.front()->nodeInfo.max.y - elementsRegions.front()->nodeInfo.min.y);
		eslog::info("  %*s : %16f %16f %16f\n", namesize - 18, "Z-COORDINATES", elementsRegions.front()->nodeInfo.min.z, elementsRegions.front()->nodeInfo.max.z, elementsRegions.front()->nodeInfo.max.z - elementsRegions.front()->nodeInfo.min.z);
		eslog::info(" ============================================================================================= \n");

		eslog::info(" %*s : %16s %16s\n", namesize, "REGION NAME", "ELEMENTS", "NODES");
		eslog::info(" ============================================================================================= \n");

		eslog::info("  %s%*s : %16s %16s\n", "TOTAL NUMBER OF NODES", namesize - 22, " ", " ", Parser::stringwithcommas(nodes->uniqInfo.totalSize).c_str());
		eslog::info("  %s%*s : %16s\n", "TOTAL NUMBER OF ELEMENTS", namesize - 25, " ", Parser::stringwithcommas(elements->distribution.process.totalSize).c_str());

		for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
			if (elements->distribution.code[etype].totalSize) {
				eslog::info(" %*s : %16s\n", namesize, ename(etype).c_str(), Parser::stringwithcommas(elements->distribution.code[etype].totalSize).c_str());
			}
		}
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n\n");

		eslog::info("  %s [%3ld]%*s : %16s %16s\n", "ELEMS REGIONS SIZES", elementsRegions.size() - 1, namesize - 26, " ", Parser::stringwithcommas(nelements).c_str(), Parser::stringwithcommas(nelementsnodes).c_str());
		eslog::info(" %*s : %16s %16s\n", namesize, elementsRegions.front()->name.c_str(), "", "");
		for (size_t r = 1; r < elementsRegions.size(); r++) {
			eslog::info(" %*s : %16s %16s\n", namesize, elementsRegions[r]->name.c_str(), Parser::stringwithcommas(elementsRegions[r]->distribution.process.totalSize).c_str(), Parser::stringwithcommas(elementsRegions[r]->nodeInfo.totalSize).c_str());
		}
		eslog::info("\n");
		eslog::info("  %s [%3d]%*s : %16s %16s\n", "FACES REGIONS SIZES", fregs, namesize - 26, " ", Parser::stringwithcommas(nfaces).c_str(), Parser::stringwithcommas(nfacenodes).c_str());
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->originalDimension == 2) {
				eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions[r]->name.c_str(), Parser::stringwithcommas(boundaryRegions[r]->distribution.process.totalSize).c_str(), Parser::stringwithcommas(boundaryRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info("\n");
		eslog::info("  %s [%3d]%*s : %16s %16s\n", "EDGES REGIONS SIZES", eregs, namesize - 26, " ", Parser::stringwithcommas(nedges).c_str(), Parser::stringwithcommas(nedgenodes).c_str());
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->originalDimension == 1) {
				eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions[r]->name.c_str(), Parser::stringwithcommas(boundaryRegions[r]->distribution.process.totalSize).c_str(), Parser::stringwithcommas(boundaryRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info("\n");
		eslog::info("  %s [%3d]%*s : %16s %16s\n", "NODES REGIONS SIZES", nregs, namesize - 26, " ", " ", Parser::stringwithcommas(nnodes).c_str());
		eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions.front()->name.c_str(), " ", " ");
		for (size_t r = 1; r < boundaryRegions.size(); r++) {
			if (boundaryRegions[r]->originalDimension == 0) {
				eslog::info(" %*s : %16s %16s\n", namesize, boundaryRegions[r]->name.c_str(), " ", Parser::stringwithcommas(boundaryRegions[r]->nodeInfo.totalSize).c_str());
			}
		}
		eslog::info(" ============================================================================================= \n");

		eslog::info("  BODY STATISTICS %22s : %16s %16s %16s\n", "ELEMENTS REGION", "ELEMENTS", "FACES", "PROPORTION");
		eslog::info(" ============================================================================================= \n");
		eslog::info("  %38s : %16s %16s %9d BODIES\n", elementsRegions[0]->name.c_str(), Parser::stringwithcommas(ecountTotal).c_str(), Parser::stringwithcommas(scountTotal).c_str(), bodies->totalSize);
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		for (size_t r = 1; r < elementsRegions.size(); r++) {
			for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
				if (elementsRegions[r]->bodyElements[b] == elementsRegions[0]->bodyElements[elementsRegions[r]->bodies[b]] && elementsRegions[r]->bodies.size() == 1) {
					eslog::info("  %38s : %16s %16s %16s\n", b ? "" : elementsRegions[r]->name.c_str(),
							Parser::stringwithcommas(elementsRegions[r]->bodyElements[b]).c_str(),
							Parser::stringwithcommas(elementsRegions[r]->bodyFaces[b]).c_str(), "BODY=REGION");
				} else {
					double ratio = (double)elementsRegions[r]->bodyElements[b] / elementsRegions[0]->bodyElements[elementsRegions[r]->bodies[b]];
					eslog::info("  %38s : %16s %16s %16f\n", b ? "" : elementsRegions[r]->name.c_str(),
							Parser::stringwithcommas(elementsRegions[r]->bodyElements[b]).c_str(),
							Parser::stringwithcommas(elementsRegions[r]->bodyFaces[b]).c_str(), ratio);
				}
			}
		}
		eslog::info(" ============================================================================================= \n");
		break;
	case OutputConfiguration::LOGGER::PARSER:
		eslog::info(" ====================================== MESH STATISTICS ====================================== \n");
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			eslog::info("mesh: region=%s, dimension=%d, elements=%d, nodes=%d\n", elementsRegions[r]->name.c_str(), dimension, elementsRegions[r]->distribution.process.totalSize, elementsRegions[r]->nodeInfo.totalSize);
		}
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			eslog::info("mesh: region=%s, dimension=%d, elements=%d, nodes=%d\n", boundaryRegions[r]->name.c_str(), boundaryRegions[r]->originalDimension, boundaryRegions[r]->distribution.process.totalSize, boundaryRegions[r]->nodeInfo.totalSize);
		}
		eslog::info("mesh: region=ALL_ELEMENTS, bodies=%d\n", bodies->totalSize);
		for (size_t r = 1; r < elementsRegions.size(); r++) {
			for (size_t b = 0; b < elementsRegions[r]->bodies.size(); ++b) {
				eslog::info("mesh: region=%s, b-elements=%d, b-faces=%d\n", elementsRegions[r]->name.c_str(), elementsRegions[r]->bodyElements[b], elementsRegions[r]->bodyFaces[b]);
			}
		}
		for (auto it = contactInterfaces.begin(); it != contactInterfaces.end(); ++it) {
			eslog::info("mesh: region=%s, s-faces=%d, s-area=%.5f, d-faces=%d, d-area=%.5f\n", (*it)->name.c_str(),
					contact->interfaces[(*it)->interfaceIndex].from.faces, contact->interfaces[(*it)->interfaceIndex].from.area,
					contact->interfaces[(*it)->interfaceIndex].to.faces, contact->interfaces[(*it)->interfaceIndex].to.area);
		}
		break;
	}
}

struct __meshinfo__ {
	struct op_min { template<typename T> T operator()(const T&i, const T&j) const { return std::min(i, j); } };
	struct op_max { template<typename T> T operator()(const T&i, const T&j) const { return std::max(i, j); } };
	struct op_sum { template<typename T> T operator()(const T&i, const T&j) const { return i + j; } };

	template <typename T>
	struct domain_t {
		T elements = 0, uniquenodes = 0, nodes = 0, neighbors = 0, ninterface = 0, einterface = 0;
		double nratio = 0, eratio = 0;

		void set(const T &value)
		{
			elements = uniquenodes = nodes = neighbors = ninterface = einterface = value;
			nratio = eratio = value;
		}

		template<typename TT>
		domain_t& operator=(const domain_t<TT> &other)
		{
			elements = other.elements;
			uniquenodes = other.uniquenodes;
			nodes = other.nodes;
			neighbors = other.neighbors;
			ninterface = other.ninterface;
			einterface = other.einterface;
			nratio = other.nratio;
			eratio = other.eratio;
			return *this;
		}

		template <typename OP>
		void apply(const domain_t &other, const OP &op)
		{
			elements = op(elements, other.elements);
			uniquenodes = op(uniquenodes, other.uniquenodes);
			nodes = op(nodes, other.nodes);
			neighbors = op(neighbors, other.neighbors);
			ninterface = op(ninterface, other.ninterface);
			einterface = op(einterface, other.einterface);
			nratio = op(nratio, other.nratio);
			eratio = op(eratio, other.eratio);
		}

		void avg(const domain_t<esint> &sum)
		{
			elements /= sum.elements;
			uniquenodes /= sum.uniquenodes;
			nodes /= sum.nodes;
			neighbors /= sum.neighbors;
			ninterface /= sum.ninterface;
			einterface /= sum.einterface;
			nratio /= sum.nratio;
			eratio /= sum.eratio;
		}

		void var(const domain_t<esint> &sum)
		{
			if (sum.elements - 1 > 0) elements = std::sqrt(elements * (1. / (sum.elements - 1)));
			if (sum.uniquenodes - 1 > 0) uniquenodes = std::sqrt(uniquenodes * (1. / (sum.uniquenodes - 1)));
			if (sum.nodes - 1 > 0) nodes = std::sqrt(nodes * (1. / (sum.nodes - 1)));
			if (sum.neighbors - 1 > 0) neighbors = std::sqrt(neighbors * (1. / (sum.neighbors - 1)));
			if (sum.ninterface - 1 > 0) ninterface = std::sqrt(ninterface * (1. / (sum.ninterface - 1)));
			if (sum.einterface - 1 > 0) einterface = std::sqrt(einterface * (1. / (sum.einterface - 1)));
			if (sum.nratio - 1 > 0) nratio = std::sqrt(nratio * (1. / (sum.nratio - 1)));
			if (sum.eratio - 1 > 0) eratio = std::sqrt(eratio * (1. / (sum.eratio - 1)));
		}
	};
	template <typename T>
	struct cluster_t: public domain_t<T> {
		T domains = 0;

		void set(const T &value)
		{
			domain_t<T>::set(value);
			domains = value;
		}

		template<typename TT>
		cluster_t& operator=(const cluster_t<TT> &other)
		{
			domain_t<T>::operator=(other);
			domains = other.domains;
			return *this;
		}

		template<typename TT>
		cluster_t& operator=(const domain_t<TT> &other)
		{
			domain_t<T>::operator=(other);
			return *this;
		}

		template <typename OP>
		void apply(const cluster_t &other, const OP &op)
		{
			domain_t<T>::apply(other, op);
			domains = op(domains, other.domains);
		}

		void avg(const cluster_t<esint> &sum)
		{
			domain_t<T>::avg(sum);
			domains /= sum.domains;
		}

		void var(const cluster_t<esint> &sum)
		{
			domain_t<T>::var(sum);
			if (sum.domains - 1 > 0) domains = std::sqrt(domains * (1. / (sum.domains - 1)));
		}
	};
	template <typename T>
	struct mpi_t: public cluster_t<T> {
		T clusters = 0;

		void set(const T &value)
		{
			cluster_t<T>::set(value);
			clusters = value;
		}

		template<typename TT>
		mpi_t& operator=(const mpi_t<TT> &other)
		{
			cluster_t<T>::operator=(other);
			clusters = other.clusters;
			return *this;
		}

		template <typename OP>
		void apply(const mpi_t &other, const OP &op)
		{
			cluster_t<T>::apply(other, op);
			clusters = op(clusters, other.clusters);
		}

		void avg(const mpi_t<esint> &sum)
		{
			cluster_t<T>::avg(sum);
			clusters /= sum.clusters;
		}

		void var(const mpi_t<esint> &sum)
		{
			cluster_t<T>::var(sum);
			if (sum.clusters - 1 > 0) clusters = std::sqrt(clusters * (1. / (sum.clusters - 1)));
		}
	};

	template<template<typename> typename C, typename T>
	struct stats_t {
		C<T> min, max, sum;
		C<double> avg;

		template<template<typename> typename CC, typename TT>
		stats_t<C, T>& operator=(const stats_t<CC, TT> &other)
		{
			min = other.min;
			max = other.max;
			sum = other.sum;
			avg = other.avg;
			return *this;
		}

		void set(const C<T> &value)
		{
			avg = min = max = sum = value;
		}

		void apply(const stats_t<C, T> &other)
		{
			min.apply(other.min, op_min());
			max.apply(other.max, op_max());
			avg.apply(other.avg, op_sum());
			sum.apply(other.sum, op_sum());
		}
	};

	struct value_t {
		mpi_t<esint> mpi;
		cluster_t<esint> cluster;
		domain_t<esint> domain;
	} value;

	struct variance_t {
		mpi_t<double> mpi;
		cluster_t<double> cluster;
		domain_t<double> domain;
	} variance;

	struct base_t {
		stats_t<mpi_t, esint> mpi;
		stats_t<cluster_t, esint> cluster;
		stats_t<domain_t, esint> domain;

		void apply(const base_t &other)
		{
			mpi.apply(other.mpi);
			cluster.apply(other.cluster);
			domain.apply(other.domain);
		}
	} stats;

	void avg()
	{
		stats.mpi.avg.avg(stats.mpi.sum);
		stats.cluster.avg.avg(stats.cluster.sum);
		stats.domain.avg.avg(stats.domain.sum);
	}

	void var()
	{
		variance.mpi.var(stats.mpi.sum);
		variance.cluster.var(stats.cluster.sum);
		variance.domain.var(stats.domain.sum);
	}

	MPI_Datatype base_mpi_t;
	MPI_Op base_mpi_op;

	static void op(void *in, void *out, int *len, MPI_Datatype *datatype)
	{
		for (int i = 0; i < *len; i++) {
			(static_cast<base_t*>(out) + i)->apply(*(static_cast<base_t*>(in) + i));
		}
	}

	__meshinfo__()
	{
		MPI_Type_contiguous(sizeof(stats), MPI_BYTE, &base_mpi_t);
		MPI_Type_commit(&base_mpi_t);
		MPI_Op_create(op, true, &base_mpi_op);
	}

	~__meshinfo__()
	{
		MPI_Type_free(&base_mpi_t);
		MPI_Op_free(&base_mpi_op);
	}
};

void Mesh::printDecompositionStatistics()
{
	__meshinfo__ mesh;
	std::vector<esint> dnodes(domains->size), dninterface(domains->size), delements(domains->size), deinterface(domains->size), dneighs(domains->size);
	std::vector<double> dnratio(domains->size), deratio(domains->size);
	std::vector<esint> cnodes(clusters->size), cninterface(clusters->size), celements(clusters->size), ceinterface(clusters->size), cneighs(clusters->size), cdomains(clusters->size);
	std::vector<double> cnratio(clusters->size), ceratio(clusters->size);

	#pragma omp parallel for
	for (int t = 0; t < info::env::threads; ++t) {
		esint inodes = 0;
		for (auto ranks = nodes->ranks->begin(t); ranks != nodes->ranks->end(t); ++ranks) {
			if (ranks->size() > 1) {
				++inodes;
			}
		}
		#pragma omp atomic
		mesh.value.mpi.ninterface += inodes;

		esint ielements = 0;
		for (auto neighs = elements->faceNeighbors->begin(t); neighs != elements->faceNeighbors->end(t); ++neighs) {
			for (auto n = neighs->begin(); n != neighs->end(); ++n) {
				if (*n != -1 && (*n < elements->distribution.process.offset || elements->distribution.process.next <= *n)) {
					++ielements;
				}
			}
		}
		#pragma omp atomic
		mesh.value.mpi.einterface += ielements;
	}

	mesh.value.mpi.clusters    = clusters->size;
	mesh.value.mpi.domains     = domains->size;
	mesh.value.mpi.elements    = elements->distribution.process.size;
	mesh.value.mpi.uniquenodes = nodes->uniqInfo.size;
	mesh.value.mpi.nodes       = nodes->size;
	mesh.value.mpi.nratio      = (double)mesh.value.mpi.ninterface / mesh.value.mpi.nodes;
	mesh.value.mpi.eratio      = (double)mesh.value.mpi.einterface / mesh.value.mpi.elements;
	mesh.value.mpi.neighbors   = neighbors.size();

	mesh.stats.mpi.set(mesh.value.mpi);
	mesh.stats.mpi.sum.set(1);

	Communication::allReduce(&mesh.stats, NULL, 1, mesh.base_mpi_t, mesh.base_mpi_op);
	mesh.avg();
	mesh.variance.mpi.elements    = std::pow(mesh.stats.mpi.avg.elements    - mesh.value.mpi.elements, 2);
	mesh.variance.mpi.uniquenodes = std::pow(mesh.stats.mpi.avg.uniquenodes - mesh.value.mpi.uniquenodes, 2);
	mesh.variance.mpi.nodes       = std::pow(mesh.stats.mpi.avg.nodes       - mesh.value.mpi.nodes, 2);
	mesh.variance.mpi.neighbors   = std::pow(mesh.stats.mpi.avg.neighbors   - mesh.value.mpi.neighbors, 2);
	mesh.variance.mpi.ninterface  = std::pow(mesh.stats.mpi.avg.ninterface  - mesh.value.mpi.ninterface, 2);
	mesh.variance.mpi.einterface  = std::pow(mesh.stats.mpi.avg.einterface  - mesh.value.mpi.einterface, 2);
	mesh.variance.mpi.nratio      = std::pow(mesh.stats.mpi.avg.nratio      - mesh.value.mpi.nratio, 2);
	mesh.variance.mpi.eratio      = std::pow(mesh.stats.mpi.avg.eratio      - mesh.value.mpi.eratio, 2);
	mesh.variance.mpi.domains     = std::pow(mesh.stats.mpi.avg.domains     - mesh.value.mpi.domains, 2);
	mesh.variance.mpi.clusters    = std::pow(mesh.stats.mpi.avg.clusters    - mesh.value.mpi.clusters, 2);
	if (mesh.stats.mpi.max.clusters == 1) {
		mesh.variance.cluster = mesh.variance.mpi;
	} else {
		for (esint c = 0; c < clusters->size; ++c) {
			mesh.variance.cluster.elements    = std::pow(mesh.stats.cluster.avg.elements    - celements[c], 2);
			mesh.variance.cluster.nodes       = std::pow(mesh.stats.cluster.avg.nodes       - cnodes[c], 2);
			mesh.variance.cluster.ninterface  = std::pow(mesh.stats.cluster.avg.ninterface  - cninterface[c], 2);
			mesh.variance.cluster.einterface  = std::pow(mesh.stats.cluster.avg.einterface  - ceinterface[c], 2);
			mesh.variance.cluster.nratio      = std::pow(mesh.stats.cluster.avg.nratio      - cnratio[c], 2);
			mesh.variance.cluster.eratio      = std::pow(mesh.stats.cluster.avg.eratio      - ceratio[c], 2);
			mesh.variance.cluster.neighbors   = std::pow(mesh.stats.cluster.avg.neighbors   - cneighs[c], 2);
		}
	}
	for (esint d = 0; d < domains->size; ++d) {
		mesh.variance.domain.elements    = std::pow(mesh.stats.domain.avg.elements    - delements[d], 2);
		mesh.variance.domain.nodes       = std::pow(mesh.stats.domain.avg.nodes       - dnodes[d], 2);
		mesh.variance.domain.ninterface  = std::pow(mesh.stats.domain.avg.ninterface  - dninterface[d], 2);
		mesh.variance.domain.einterface  = std::pow(mesh.stats.domain.avg.einterface  - deinterface[d], 2);
		mesh.variance.domain.nratio      = std::pow(mesh.stats.domain.avg.nratio      - dnratio[d], 2);
		mesh.variance.domain.eratio      = std::pow(mesh.stats.domain.avg.eratio      - deratio[d], 2);
		mesh.variance.domain.neighbors   = std::pow(mesh.stats.domain.avg.neighbors   - dneighs[d], 2);
	}

	Communication::allReduce(&mesh.variance, NULL, sizeof(mesh.variance) / sizeof(double), MPI_DOUBLE, MPI_SUM);
	mesh.var();

	eslog::info("\n ================================== DECOMPOSITION STATISTICS ================================= \n");
	if (info::config::output.logger == OutputConfiguration::LOGGER::PARSER) {
		eslog::info("decomposition: clusters per MPI: %d\n", mesh.stats.mpi.max.clusters);
		return;
	}
	if (mesh.stats.mpi.max.clusters > 1) {
		eslog::info(" ============================================================================================= \n");
		eslog::info("  PER MPI STATISTICS         %12s %12s %12s %12s %12s \n", "MIN", "MAX", "AVG", "VAR", "IMBALANCE");
		eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
		eslog::info("  CLUSTERS                 : %12d %12d %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.clusters,
				mesh.stats.mpi.max.clusters,
				mesh.stats.mpi.avg.clusters,
				mesh.variance.mpi.clusters,
				mesh.stats.mpi.min.clusters ? (double)mesh.stats.mpi.max.clusters / mesh.stats.mpi.min.clusters : 0.0);
		eslog::info("  DOMAINS                  : %12d %12d %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.domains,
				mesh.stats.mpi.max.domains,
				mesh.stats.mpi.avg.domains,
				mesh.variance.mpi.domains,
				mesh.stats.mpi.min.domains ? (double)mesh.stats.mpi.max.domains / mesh.stats.mpi.min.domains : 0.0);
		eslog::info("  UNIQUE NODES             : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.uniquenodes).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.uniquenodes).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.uniquenodes).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.uniquenodes).c_str(),
				mesh.stats.mpi.min.uniquenodes ? (double)mesh.stats.mpi.max.uniquenodes / mesh.stats.mpi.min.uniquenodes : 0.0);
		eslog::info("  NODES                    : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.nodes).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.nodes).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.nodes).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.nodes).c_str(),
				mesh.stats.mpi.min.nodes ? (double)mesh.stats.mpi.max.nodes / mesh.stats.mpi.min.nodes : 0.0);
		eslog::info("  INTERFACE NODES          : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.ninterface).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.ninterface).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.ninterface).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.ninterface).c_str(),
				mesh.stats.mpi.min.ninterface ? (double)mesh.stats.mpi.max.ninterface / mesh.stats.mpi.min.ninterface : 0.0);
		eslog::info("  INTERFACE NODES RATIO    : %12.2f %12.2f %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.nratio,
				mesh.stats.mpi.max.nratio,
				mesh.stats.mpi.avg.nratio,
				mesh.variance.mpi.nratio,
				mesh.stats.mpi.min.nratio ? mesh.stats.mpi.max.nratio / mesh.stats.mpi.min.nratio : 0.0);
		eslog::info("  ELEMENTS                 : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.elements).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.elements).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.elements).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.elements).c_str(),
				mesh.stats.mpi.min.elements ? (double)mesh.stats.mpi.max.elements / mesh.stats.mpi.min.elements : 0.0);
		eslog::info("  INTERFACE ELEMENTS       : %12s %12s %12s %12s %12.2f \n",
				Parser::stringwithcommas(mesh.stats.mpi.min.einterface).c_str(),
				Parser::stringwithcommas(mesh.stats.mpi.max.einterface).c_str(),
				Parser::stringwithcommas((int)mesh.stats.mpi.avg.einterface).c_str(),
				Parser::stringwithcommas((int)mesh.variance.mpi.einterface).c_str(),
				mesh.stats.mpi.min.einterface ? (double)mesh.stats.mpi.max.einterface / mesh.stats.mpi.min.einterface : 0);
		eslog::info("  INTERFACE ELEMENTS RATIO : %12.2f %12.2f %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.eratio,
				mesh.stats.mpi.max.eratio,
				mesh.stats.mpi.avg.eratio,
				mesh.variance.mpi.eratio,
				mesh.stats.mpi.min.eratio ? mesh.stats.mpi.max.eratio / mesh.stats.mpi.min.eratio : 0.0);
		eslog::info("  NEIGHBORS                : %12d %12d %12.2f %12.2f %12.2f \n",
				mesh.stats.mpi.min.neighbors,
				mesh.stats.mpi.max.neighbors,
				mesh.stats.mpi.avg.neighbors,
				mesh.variance.mpi.neighbors,
				mesh.stats.mpi.min.neighbors ? (double)mesh.stats.mpi.max.neighbors / mesh.stats.mpi.min.neighbors : 0.0);
	}
	eslog::info(" ============================================================================================= \n");
//
//	if (info::config::output.store_decomposition > 1) {
//		NamedData *data;
//		if (mesh.stats.mpi.max.clusters > 1) {
//			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_CLUSTERS");
//			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.clusters);
//			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_NODES");
//			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.nodes);
//			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_INODES");
//			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.ninterface);
//			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_INODES_RATIO");
//			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.nratio);
//			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_ELEMENTS");
//			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.elements);
//			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_IELEMENTS");
//			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.einterface);
//			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_IELEMENTS_RATIO");
//			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.eratio);
//			data = info::mesh->elements->appendData(1, NamedData::DataType::SCALAR, "MPI_NEIGHBORS");
//			std::fill(data->data.begin(), data->data.end(), mesh.value.mpi.neighbors);
//		}
//	}
}

