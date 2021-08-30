
#include "insitu.h"

namespace mesio {

struct Catalyst {};

InSitu::InSitu()
: _catalyst(NULL)
{

}

InSitu::~InSitu()
{
	if (_catalyst != NULL) {
		delete _catalyst;
	}
}

void InSitu::updateMesh()
{
	_catalyst = new Catalyst();
}

void InSitu::updateSolution()
{
//	_catalyst->update();
//	sleep(info::config::output.catalyst_sleep_time);
}

}


