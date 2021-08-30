
#include "eslog.hpp"
#include "basis/logging/logger.h"
#include "basis/logging/timelogger.h"
#include "wrappers/mpi/communication.h"
#include "basis/utilities/sysutils.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/systeminfo.h"
#include "esinfo/config.h"
#include "esinfo/envinfo.h"

namespace mesio {
namespace eslog {

char buffer[BUFFER_SIZE];
LoggerBase *logger = NULL;

double time()
{
	return TimeLogger::time();
}

double duration()
{
	return TimeLogger::duration();
}

void init(LoggerBase *logger)
{
	eslog::logger = logger;

	for (int i = 0; i < logger->size; i++) {
		logger->args[i]->rank = info::mpi::rank;
		logger->args[i]->grank = info::mpi::grank;
	}
	Communication::broadcast(&TimeLogger::initTime, sizeof(time_t), MPI_BYTE, 0, MPITools::global);
}

void initFiles()
{
	logger->initOutput();
}

void reinit()
{
	for (int i = 0; i < logger->size; i++) {
		logger->args[i]->rank = info::mpi::rank;
		logger->args[i]->grank = info::mpi::grank;
	}
}

void printRunInfo(int *argc, char ***argv)
{
	int width = 78;

	struct tm *timeinfo;
	char date[80], time[80];
	timeinfo = std::localtime(&TimeLogger::initTime);
	std::strftime(date, 80, "%F", timeinfo);
	std::strftime(time, 80, "%H-%M-%S", timeinfo);

	auto divide = [] (std::string &str, size_t max, const char* separator = "/") {
		std::string backup(str);
		size_t i = 0, c = 0, cmax = 10;
		while (++c < cmax && str.size() - i > max) {
			size_t it = str.find_last_of(separator, i + max);
			str.insert(it, " ==\n == ");
			size_t spaces = i + max - it + (i ? 9 : 0);
			str.insert(i, spaces, ' ');
			i = it + 8 + spaces;
		}
		str.insert(i, max - (str.size() - i) + (i ? 9 : 0), ' ');
		if (c == cmax) {
			str = backup;
		}
	};

	std::string cxx = info::system::cxx();
	std::string mesh = info::config::input.path;
	std::string runpath(info::env::pwd() ? info::env::pwd() : "");
	std::string cmd((*argv)[0]);
	if (runpath.size()) {
		if (mesh[0] != '/') {
			mesh = (runpath + "/" + mesh);
		}
	}

	divide(cxx, width);
	divide(mesh, width);
	divide(runpath, width);
	divide(cmd, width, ":/ ");

	eslog::info(" ======================================= MESIO RUN INFO ====================================== \n");
	eslog::info(" ============================================================================================= \n");
	eslog::info(" == CXX      %*s == \n", width, cxx.c_str());
	eslog::info(" == CXXFLAGS %*s == \n", width, info::system::cxxflags());
	eslog::info(" == COMMIT   %*s == \n", width, info::system::commit());
	eslog::info(" == DATE [YYYY-MM-DD]   %*s == \n", width - 11, date);
	eslog::info(" == TIME [HH-MM-SS]    %*s == \n", width - 10, time);
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
	eslog::info(" == COMMAND  %*s == \n", width, cmd.c_str());
	eslog::info(" == MPI_COMM_WORLD %*d == \n", width - 6, info::mpi::size);
	eslog::info(" == OMP_NUM_THREADS %*d == \n", width - 7, info::env::threads);
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
	eslog::info(" == NUMBER OF LOADERS %*d == \n", width - 9, MPITools::subset->acrosssize);
	eslog::info(" == NUMBER OF WRITTERS %*d == \n", width - 10, MPITools::subset->acrosssize);
	switch (info::config::output.mode) {
	case OutputConfiguration::MODE::SYNC  : eslog::info(" == STORING MODE %*s == \n", width - 4, "SYNCHRONIZED"); break;
	case OutputConfiguration::MODE::PTHREAD: eslog::info(" == STORING MODE %*s == \n", width - 4, "SEPARETED P-THREAD"); break;
	}
	eslog::info(" ==    -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -    == \n");
	eslog::info(" == RUN PATH %*s == \n", width, runpath.c_str());
	eslog::info(" ============================================================================================= \n");
}

void finish()
{
	eslog::info("\n ======================================== RUN FINISHED ======================================= \n");
	logger->finish();
	if (logger) delete logger;
}

void always()
{
	logger->always();
}

void start(const char* name, const char* section)
{
	logger->start(name, section);
}

void checkpoint(const char* name)
{
	logger->checkpoint(name);
}

void end(const char* name)
{
	logger->end(name);
}

void ln()
{
	logger->ln();
}

void nextStep(int step)
{
	logger->nextLoadStep(step);
	eslog::info("\n ============================================================================================= \n");
	eslog::info(" ======================================= LOAD STEP %3d ======================================= \n", step);
	eslog::info(" ============================================================================================= \n");
}

void startln(const char* name, const char* section)
{
	start(name, section);
	ln();
}

void checkpointln(const char* name)
{
	checkpoint(name);
	ln();
}

void endln(const char* name)
{
	end(name);
	logger->ln();
}

void param(const char* name, const int &value)
{
	logger->param(name, value);
}

void param(const char* name, const long &value)
{
	logger->param(name, value);
}

void param(const char* name, const long unsigned int &value)
{
	logger->param(name, value);
}

void param(const char* name, const double &value)
{
	logger->param(name, value);
}

void param(const char* name, const char* value)
{
	logger->param(name, value);
}

void info(const char* msg)
{
	logger->output(msg, VerboseArg::COLOR::WHITE);
}

void solver(const char* msg)
{
	logger->output(msg, VerboseArg::COLOR::WHITE);
}

void linearsolver(const char* msg)
{
	for (int n = 0; n < logger->size; ++n) {
		if (logger->args[n]->argflag == 'v' && logger->args[n]->verbosity > 2) {
			logger->output(msg, VerboseArg::COLOR::WHITE);
			break;
		}
	}
}

void duration(const char* msg)
{
	for (int n = 0; n < logger->size; ++n) {
		if (logger->args[n]->argflag == 'm' && logger->args[n]->verbosity > 1) {
			logger->output(msg, VerboseArg::COLOR::WHITE);
			break;
		}
	}
}

void warning(const char* msg)
{
	logger->output(msg, VerboseArg::COLOR::YELLOW);
}

void storedata(const char* msg)
{
	logger->output(msg, VerboseArg::COLOR::BLUE);
}

void failure(const char* msg)
{
	logger->error(msg);
	utils::printStack();
	fflush(stderr);
	exit(EXIT_FAILURE);
}

void internalFailure(const char* msg)
{
	logger->error("INTERNAL FAILURE: ");
	logger->error(msg);
	utils::printStack();
	fflush(stderr);
	exit(EXIT_FAILURE);
}

void error(const char* msg)
{
	logger->error(msg);
	fflush(stderr);
	exit(EXIT_FAILURE);
}

void globalerror(const char* msg)
{
	if (info::mpi::rank == 0) {
		logger->error(msg);
	}
	Communication::barrier();
	MPI_Finalize();
	exit(EXIT_FAILURE);
}

}
}



