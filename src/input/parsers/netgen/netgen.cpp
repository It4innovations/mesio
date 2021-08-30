
#include "netgen.h"
#include "parser/neutralmesh.h"
#include "config/input.h"
#include "esinfo/eslog.h"
#include "input/meshbuilder.h"
#include "input/parsers/distributedscanner.h"

using namespace mesio;

NetgenNeutralLoader::NetgenNeutralLoader(const InputConfiguration &configuration)
: _configuration(configuration)
{


}

void NetgenNeutralLoader::load()
{
	eslog::startln("NETGEN PARSER: STARTED", "NETGEN PARSER");

	InputFilePack meshfile;
	meshfile.commitFiles({ _configuration.path });
	meshfile.prepare();
	eslog::checkpointln("NETGEN PARSER: MESH READER PREPARED");

	meshfile.read();
	eslog::checkpointln("NETGEN PARSER: MESH READ");

	meshfile.next();
	DistributedScanner::align(meshfile, "\n");

	NetgenNeutralMesh mesh(meshfile);

	mesh.parse(*this);
	body.resize(etype.size());
	material.resize(etype.size());
	eslog::endln("NETGEN PARSER: GEOMETRY PARSED");
}
