
#include "xdmf.h"
#include "lightdata/lightdata.h"
#include "heavydata/griddata.h"
#include "esinfo/eslog.h"
#include "config/input.h"

using namespace mesio;


XDMFLoader::XDMFLoader(const InputConfiguration &configuration)
: _configuration(configuration)
{

}

void XDMFLoader::load()
{
	eslog::startln("XDMF PARSER: STARTED", "XDMF PARSER");

	LightData lightdata(_configuration.path);
	eslog::checkpointln("XDMF PARSER: LIGHTDATA PARSED");

	GridData grid(lightdata);
	grid.scan();
	eslog::checkpointln("XDMF PARSER: LIGHTDATA SCANNED");

	grid.read();
	eslog::checkpointln("XDMF PARSER: HEAVYDATA READ");

	grid.parse(*this);
	removeDuplicates = true;
	body.resize(etype.size());
	material.resize(etype.size());
	eslog::endln("XDMF PARSER: HEAVYDATA PARSED");
}
