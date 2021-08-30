
#include "ensight.h"
#include "parser/casefile.h"
#include "parser/geometry.h"

#include "config/input.h"
#include "esinfo/eslog.h"
#include "esinfo/mpiinfo.h"
#include "input/meshbuilder.h"

using namespace mesio;

EnsightLoader::EnsightLoader(const InputConfiguration &configuration)
: _configuration(configuration)
{

}

void EnsightLoader::load()
{
	eslog::startln("ENSIGHT PARSER: STARTED", "ENSIGHT PARSER");

	EnsightCasefile casefile(_configuration.path);
	eslog::checkpointln("ENSIGHT PARSER: CASEFILE READ");

	InputFilePack geofile;
	geofile.commitFiles({ casefile.geometry });
	geofile.prepare();
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY READER PREPARED");

	geofile.read();
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY READ");

	geofile.next();

	EnsightGeometry geometry(geofile);
	geometry.scan();
	eslog::checkpointln("ENSIGHT PARSER: GEOMETRY SCANNED");

	geometry.parse(*this);
	removeDuplicates = true;
	body.resize(etype.size());
	material.resize(etype.size());
	eslog::endln("ENSIGHT PARSER: GEOMETRY PARSED");
}
