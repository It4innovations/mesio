
#ifndef SRC_CONFIG_ECF_INPUT_INPUT_H_
#define SRC_CONFIG_ECF_INPUT_INPUT_H_

#include <string>
#include <map>

namespace mesio {

struct ParMETISConfiguration {
	bool refinement = false;
	double tolerance = 1.05;
};

struct METISConfiguration {

	enum class OBJECTIVE_TYPE {
		VOLUME, EDGECUT
	};

	OBJECTIVE_TYPE objective_type = OBJECTIVE_TYPE::VOLUME;
	int continuous = 1;
};

struct PTScotchConfiguration {

};

struct ScotchConfiguration {

};

struct KaHIPConfiguration {

};

struct DecompositionConfiguration {

	enum class ParallelDecomposer {
		NONE,
		METIS,
		PARMETIS,
		PTSCOTCH,
		HILBERT_CURVE
	};

	enum class SequentialDecomposer {
		NONE,
		METIS,
		SCOTCH,
		KAHIP
	};

	ParallelDecomposer parallel_decomposer = ParallelDecomposer::NONE;
	SequentialDecomposer sequential_decomposer = SequentialDecomposer::NONE;
	int mesh_duplication = 1;
	int domains = 0;

	bool force_continuity = 0;
	bool separate_materials = 0, separate_regions = 0, separate_etypes = 0;
	ParMETISConfiguration parmetis_options;
	METISConfiguration metis_options;
	PTScotchConfiguration ptscotch_options;
	ScotchConfiguration scotch_options;
	KaHIPConfiguration kahip_options;
};

struct InputTransformationConfiguration {

	enum class TRANSFORMATION {
		TRANSLATE,
		ROTATE,
		SCALE,
		SHEAR
	};

	TRANSFORMATION transformation;
	double x, y, z;

	int instances;
};

struct InputConfiguration {

	enum class FORMAT {
		ANSYS_CDB,
		OPENFOAM,
		ABAQUS,
		XDMF,
		ENSIGHT,
		VTK_LEGACY,
		NETGET,
		NEPER
	};

	enum class LOADER {
		MPI,
		MPI_COLLECTIVE,
		POSIX
	};

	std::string path;
	FORMAT format = FORMAT::ANSYS_CDB;

	bool omit_midpoints = false, insert_midpoints = false;
	bool omit_face_sets = false;
	bool keep_material_sets = false;
//	bool convert_database;
	double duplication_tolerance = 1e-6;

//	bool insert_orientation;

	LOADER loader = LOADER::POSIX;
	size_t stripe_size = 1024 * 1024;
	int third_party_scalability_limit = 1024;

	DecompositionConfiguration decomposition;
	std::map<std::string, InputTransformationConfiguration> transformations;
};

}



#endif /* SRC_CONFIG_ECF_INPUT_INPUT_H_ */
