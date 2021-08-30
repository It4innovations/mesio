
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_

#include "visualization.h"

namespace mesio {

class Catalyst;

struct InSitu: public Visualization {

	InSitu();
	~InSitu();

	void updateMesh();
	void updateSolution();

protected:
	Catalyst *_catalyst;
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_SEPARATED_INSITU_H_ */
