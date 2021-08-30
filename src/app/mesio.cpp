
#include "esinfo/config.h"
#include "esinfo/eslog.hpp"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/systeminfo.h"
#include "wrappers/mpi/communication.h"

#include "basis/logging/logger.h"
#include "basis/logging/progresslogger.h"
#include "basis/logging/timelogger.h"

#include "output/output.h"

#include <getopt.h>
#include <cstring>

using namespace mesio;

bool set(int &argc, char** &argv)
{
	int c, set = 0;
	while ((c = getopt (argc, argv, "p:i:o:s:")) != -1)
		switch (c) {
		case 'p':
			set |= 1;
			info::config::input.path = optarg;
			break;
		case 'i':
			if (memcmp(optarg, "ANSYS_CDB", 9) == 0) {
				set |= 2;
				info::config::input.format = InputConfiguration::FORMAT::ANSYS_CDB;
			}
			if (memcmp(optarg, "OPENFOAM", 8) == 0) {
				set |= 2;
				info::config::input.format = InputConfiguration::FORMAT::OPENFOAM;
			}
			if (memcmp(optarg, "ABAQUS", 6) == 0) {
				set |= 2;
				info::config::input.format = InputConfiguration::FORMAT::ABAQUS;
			}
			if (memcmp(optarg, "XDMF", 4) == 0) {
				set |= 2;
				info::config::input.format = InputConfiguration::FORMAT::XDMF;
			}
			if (memcmp(optarg, "ENSIGHT", 7) == 0) {
				set |= 2;
				info::config::input.format = InputConfiguration::FORMAT::ENSIGHT;
			}
			if (memcmp(optarg, "VTK_LEGACY", 10) == 0) {
				set |= 2;
				info::config::input.format = InputConfiguration::FORMAT::VTK_LEGACY;
			}
			if (memcmp(optarg, "NETGET", 6) == 0) {
				set |= 2;
				info::config::input.format = InputConfiguration::FORMAT::NETGET;
			}
			if (memcmp(optarg, "NEPER", 5) == 0) {
				set |= 2;
				info::config::input.format = InputConfiguration::FORMAT::NEPER;
			}
			break;
		case 'o':
			if (memcmp(optarg, "VTK_LEGACY", 10) == 0) {
				set |= 4;
				info::config::output.format = OutputConfiguration::FORMAT::VTK_LEGACY;
			}
			if (memcmp(optarg, "ENSIGHT", 7) == 0) {
				set |= 4;
				info::config::output.format = OutputConfiguration::FORMAT::ENSIGHT;
			}
			if (memcmp(optarg, "XDMF", 4) == 0) {
				set |= 4;
				info::config::output.format = OutputConfiguration::FORMAT::XDMF;
			}
			if (memcmp(optarg, "STL_SURFACE", 11) == 0) {
				set |= 4;
				info::config::output.format = OutputConfiguration::FORMAT::STL_SURFACE;
			}
			if (memcmp(optarg, "NETGEN", 6) == 0) {
				set |= 4;
				info::config::output.format = OutputConfiguration::FORMAT::NETGEN;
			}
			break;
		case 's':
			set |= 8;
			info::config::output.path = optarg;
			break;
		case '?':
			if (optopt == 'p' || optopt == 'i' || optopt == 'o' || optopt == 's') {
				eslog::info(" MESIO: Option -%c requires an argument.\n", optopt);
			} else {
				eslog::info(" MESIO: Unknown option character `\\x%x'.\n", optopt);
			}
			break;
		default:
			break;
	}

	if (set != 15) {
		eslog::info(" MESIO: INVALID CONFIGURATION DETECTED.\n");
		eslog::info(" MESIO:\n");
		eslog::info(" MESIO: RUN:\n");
		eslog::info(" MESIO: mpirun -n #PROCS mesio -i INPUT_FORMAT -p INPUT_PATH -o OUTPUT_FORMAT -s STORE_PATH\n");
		eslog::info(" MESIO:\n");
		eslog::info(" MESIO: FINISHED\n");
		return false;
	}
	return true;
}

int main(int argc, char** argv)
{
	info::system::setSignals();
	info::env::set();
	info::mpi::init(&argc, &argv);
	MPITools::init();

	eslog::init(new Logger<TimeLogger, ProgressTerminalLogger>);
	eslog::startln("MESIO: STARTED", "MESIO");

	if (set(argc, argv)) {
//		MPITools::setSubset(info::config::input.third_party_scalability_limit);
		info::config::output.mode = OutputConfiguration::MODE::SYNC;

		eslog::printRunInfo(&argc, &argv);
		Mesh::init();
		eslog::checkpointln("MESIO: RUN INITIALIZED");

		Mesh::load();
		eslog::checkpointln("MESIO: MESH LOADED");

		info::mesh->preprocess();
		eslog::checkpointln("MESIO: MESH PREPROCESSED");
	//	info::mesh->printMeshStatistics();
	//	info::mesh->printDecompositionStatistics();

		info::mesh->output->updateMesh();
		eslog::endln("MESIO: MESH STORED");

		Mesh::finish();
		eslog::finish();
	}
	MPITools::finish();

	info::mpi::finish();

	return 0;
}
