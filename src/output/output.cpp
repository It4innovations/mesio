
#include "output.h"

#include "visualization/visualization.h"
#include "visualization/vtklegacy.h"
#include "visualization/ensightgold.h"
#include "visualization/xdmf.h"
#include "visualization/stl.h"
#include "visualization/netgen.h"
#include "visualization/insitu.h"

#include "basis/utilities/sysutils.h"
#include "esinfo/config.h"
#include "esinfo/eslog.h"
#include "esinfo/meshinfo.h"
#include "wrappers/pthread/w.pthread.h"

#include <vector>

namespace mesio {

struct OutputExecutor {
	void insert(OutputWriter *writer)
	{
		writers.push_back(writer);
	}

	virtual ~OutputExecutor()
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			delete writers[i];
		}
	}

	virtual void mesh() = 0;
	virtual void solution() = 0;

	std::vector<OutputWriter*> writers;
};

class DirectOutputExecutor: public OutputExecutor {
public:
	virtual void mesh()
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			writers[i]->updateMesh();
		}
	}

	virtual void solution()
	{
		for (size_t i = 0; i < writers.size(); ++i) {
			writers[i]->updateSolution();
		}
	}
};



class AsyncOutputExecutor: public DirectOutputExecutor, public Pthread::Executor, public Pthread {
	struct SharedData {
		enum class TAG {
			MESH,
			SOLUTION,
			DUMMY
		} tag;

		SharedData(TAG tag): tag(tag) {}
	} app, thread;

public:
	AsyncOutputExecutor(): Pthread(this), app(SharedData::TAG::DUMMY), thread(app)
	{

	}

	void copy()
	{
		info::mesh->toBuffer();
		thread = app;
	}

	void call()
	{
		if (thread.tag == SharedData::TAG::MESH) {
			DirectOutputExecutor::mesh();
		}
		if (thread.tag == SharedData::TAG::SOLUTION) {
			DirectOutputExecutor::solution();
		}
	}

	virtual void mesh()
	{
		app = SharedData(SharedData::TAG::MESH); // dummy type
		Pthread::call();
	}

	virtual void solution()
	{
		app = SharedData(SharedData::TAG::SOLUTION);
		Pthread::call();
	}
};

}

using namespace mesio;

OutputWriter::OutputWriter()
: _path(info::config::output.path + "/"), _directory("PREPOSTDATA/"),
  _measure(info::config::output.mode == OutputConfiguration::MODE::SYNC), _allowed(true)
{
	size_t namebegin = info::config::input.path.find_last_of("/") + 1;
	size_t nameend = info::config::input.path.find_last_of(".");
	_name = info::config::input.path.substr(namebegin, nameend - namebegin);
	createOutputDirectory();
}

void OutputWriter::createOutputDirectory()
{
	utils::createDirectory({ _path, _directory });
}

Output::Output()
: _direct(new DirectOutputExecutor()), _async(NULL)
{
	if (info::config::output.mode != OutputConfiguration::MODE::SYNC) {
		_async = new AsyncOutputExecutor();
	}
	if (true) {
		OutputWriter *writer = NULL;
		switch (info::config::output.format) {
		case OutputConfiguration::FORMAT::VTK_LEGACY: writer = new VTKLegacy(); break;
		case OutputConfiguration::FORMAT::ENSIGHT: writer = new EnSightGold(); break;
		case OutputConfiguration::FORMAT::XDMF: writer = new XDMF(); break;
		case OutputConfiguration::FORMAT::STL_SURFACE: writer = new STL(); break;
		case OutputConfiguration::FORMAT::NETGEN: writer = new Netgen(); break;
		default:
			eslog::internalFailure("implement the selected output format.\n");
		}
		switch (info::config::output.mode) {
		case OutputConfiguration::MODE::SYNC: _direct->insert(writer); break;
		case OutputConfiguration::MODE::PTHREAD: _async->insert(writer); break;
		default:
			eslog::internalFailure("implement the selected output mode.\n");
		}
	}

	if (_direct->writers.size() == 0) {
		delete _direct;
		_direct = nullptr;
	}

	if (_async && _async->writers.size() == 0) {
		delete _async;
		_async = nullptr;
	}
}

Output::~Output()
{
	if (_direct) { delete _direct; }
	if (_async) { delete _async; }
}

void Output::updateMesh()
{
	// Time is not used during storing the solution
	if (_async) { _async->mesh(); }
	if (_direct) { _direct->mesh(); }
}

void Output::updateSolution()
{
	if (_allowed) {
		if (_async) { _async->solution(); }
		if (_direct) { _direct->solution(); }
	}
}

void Output::suppress()
{
	_allowed = false;
}

void Output::permit()
{
	_allowed = true;
}
