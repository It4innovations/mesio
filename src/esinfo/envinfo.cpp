
#include "envinfo.h"

#include <sstream>
#include <cstdlib>
#include <omp.h>

int mesio::info::env::threads = 1;

void mesio::info::env::set()
{
	auto getEnv = [] (int &value, const char *name)
	{
		char *var = getenv(name);
		if (var != NULL) {
			std::stringstream ss(var);
			ss >> value;
		}
	};

	getEnv(threads, "threads");

	omp_set_num_threads(threads);
}

char* mesio::info::env::pwd()
{
	return getenv("PWD");
}


