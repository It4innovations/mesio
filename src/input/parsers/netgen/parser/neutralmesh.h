
#ifndef SRC_INPUT_PARSERS_NETGEN_PARSER_NEUTRALMESH_H_
#define SRC_INPUT_PARSERS_NETGEN_PARSER_NEUTRALMESH_H_

#include "basis/io/inputfile.h"

namespace mesio {

struct MeshBuilder;

class NetgenNeutralMesh {

public:
	NetgenNeutralMesh(InputFilePack &meshfile);

	void parse(MeshBuilder &mesh);

protected:
	InputFilePack &_meshfile;
};
}



#endif /* SRC_INPUT_PARSERS_NETGEN_PARSER_NEUTRALMESH_H_ */
