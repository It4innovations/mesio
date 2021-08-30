
#include "visualization.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/config.h"
#include "mesh/store/nameddata.h"

using namespace mesio;

Visualization::Visualization()
{
//	if (info::config::output.results_store_frequency != OutputConfiguration::STORE_FREQUENCY::NEVER) {
////		createOutputDirectory();
//	}
}

Visualization::~Visualization()
{

}

bool Visualization::isRoot()
{
	return info::mpi::rank == 0;
//	if (info::mpi::rank == 0) {
//		if (step::type == step::TYPE::FTT) {
//			return true;
//		}
//		if (step::duplicate::instances == 1 && info::mpi::grank == 0) {
//			return true;
//		}
//	}
//	return false;
}

bool Visualization::storeStep()
{
	return false;
}

bool Visualization::storeData(const NamedData *data)
{
	return false;
}

Point Visualization::shrink(const Point &p, const Point &ccenter, const Point &dcenter, double cratio, double dratio) {
	Point point = ccenter + (p - ccenter) * cratio;
	point = dcenter + (point - dcenter) * dratio;
	return point;
}
