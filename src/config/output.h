
#ifndef SRC_CONFIG_ECF_OUTPUT_H_
#define SRC_CONFIG_ECF_OUTPUT_H_

namespace mesio {

struct OutputConfiguration {

	enum class FORMAT {
		VTK_LEGACY = 0,
		ENSIGHT,
		XDMF,
		STL_SURFACE,
		NETGEN
	};

	enum class WRITER {
		MPI,
		MPI_COLLECTIVE,
		POSIX
	};

	enum class LOGGER {
		USER,
		PARSER
	};

	enum class MODE {
		SYNC,
		PTHREAD,
	};

	size_t verbose_level = 1;
	size_t measure_level = 0;
	LOGGER logger = LOGGER::USER;

	FORMAT format = FORMAT::ENSIGHT;
	MODE mode = MODE::SYNC;

	WRITER writer = WRITER::MPI_COLLECTIVE;
	size_t stripe_size = 1024 * 1024, stripe_count = 1;
	int debug = 0;

	std::string path = ".";
};

}



#endif /* SRC_CONFIG_ECF_OUTPUT_H_ */
