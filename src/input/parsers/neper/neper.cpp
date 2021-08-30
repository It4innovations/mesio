
#include "neper.h"
#include "parser/msh.parser.h"
#include "config/input.h"
#include "esinfo/eslog.h"
#include "input/meshbuilder.h"
#include "input/parsers/distributedscanner.h"

using namespace mesio;

NeperLoader::NeperLoader(InputConfiguration &configuration)
: _configuration(configuration)
{
//	_configuration.insert_orientation = true;
}

void NeperLoader::load()
{
	eslog::startln("NEPER PARSER: STARTED", "NEPER PARSER");

	InputFilePack meshfile;
	meshfile.commitFiles({ _configuration.path });
	meshfile.prepare();
	eslog::checkpointln("NEPER PARSER: MESH READER PREPARED");

	meshfile.read();
	eslog::checkpointln("NEPER PARSER: MESH READ");

	meshfile.next();
	DistributedScanner::align(meshfile, "\n");

	NeperMshMesh mesh(meshfile);

	mesh.parse(*this);
	body.resize(etype.size());
	material.resize(etype.size());
	eslog::endln("NEPER PARSER: GEOMETRY PARSED");
}
