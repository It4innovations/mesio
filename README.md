
Mesio is a highly-parallel loader and converter of external unstructured meshes databases. It is composed of several commonly used algorithms that together provide a robust solution for the pre/post-processing phase of solving large-scale engineering problems ([pdf](https://doi.org/10.1109/IPDPS.2019.00084)). With mesio a user is able to use the same sequential database file with an arbitrary number of MPI processes. It allows engineers to use their favorite tools for creation of numerical models without any penalty for a parallel run.

Mesio is composed from (i) the mesh builder that is able to reconstruct an unstructured mesh from randomly scattered data across parallel processes in a distributed memory of an HPC machine without gathering data into a single node and (ii) lightweight parallel parsers of given mesh databases. Hence, the usage is not limited to any particular format -- only a simple parallel parser has to be implemented for a new format. Currently, the following parsers are implemented:
 - ANSYS CDB
 - Ensight
 - VTK Legacy
 - XDMF
 - Netgen
 - Neper
 - OpenFOAM (partially)
 - Abaqus (partially)

An output database stored by mesio is also in a sequential form for simple by a favorite visualization tool. The following format are available:
 - VTK Legacy
 - XDMF
 - Ensight
 - STL surface

Mesio functionality is provided to other researchers by [API](#mesio-api).

---
---
---

## Installation
---

####  External Dependencies

In the current version the modules can be compiled and executed on Linux operating systems only. Some functionality requires third party libraries that should be installed in the system before the installation process. Currently available wrappers are the followings:
1. Parallel graph partitioners:
 - [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
 - [PT-Scotch](https://www.labri.fr/perso/pelegrin/scotch/)
2. Sequential graph partitioners:
 - [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
 - [Scotch](https://www.labri.fr/perso/pelegrin/scotch/)
 - [KaHIP](http://algo2.iti.kit.edu/kahip/)
3. Other libraries:
 - [HDF5](https://www.hdfgroup.org/solutions/hdf5/)

#### Building the library

For compilation Python-based framework [Waf](https://waf.io/book/) is used. The compilation process has two phases: **configuration** and **compilation**.  The former configures persistent data and checks available headers and libraries. The latter builds the library. It is mandatory to re-run the configuration process after any environment change. The following commands build all modules if require libraries are available:
```sh
$ ./waf configure
$ ./waf
```
The compilation process builds all libraries and executable tools into the *build* directory. This directory should be added to ``LD_LIBRARY_PATH`` and ``PATH`` environment variables. Then it is possible to run `mesio` by the following command:
```sh
$ mpirun -n $N mesio -i INPUT_FORMAT -p INPUT_PATH -o OUTPUT_FORMAT -s STORE_PATH
```
where **N** is the number of MPI processes.

#### Set up the environment

Before running the library, the following variables should be set according to the number of CPU cores per compute node (nCores) and number of MPI processes processed per node (PPN):

 - OMP_NUM_THREADS - should be set to nCores/PPN


---
---
---

## MESIO API

Our parallel loader can be utilized by third party software by provided C API. The interface is available via the `mesio.h` header in the include directory and the `libmesioapi` library. An example of calling the library can be found in `src/api/api.mesio.cpp`. The code below shows how the loader should be called. Once the method `MESIOLoad` is finished an input database is loaded. Then, one can use provided functions to return mesh data stored in internal structures. The code shows requesting of nodes and elements only. For the full API description see the provided example and our wiki (the wiki also contains the description of how to implement a simple parallel parser).

```cpp
#include "mesio.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	MESIO mesio;
	int verbose = 2;
	int domains = 0;
	MESIOInit(MPI_COMM_WORLD, verbose);
	MESIOLoad(&mesio, MESIO_ANSYS, "path_to_file", MESIO_PARMETIS, domains);

	{ // GET NODES
		MESIOInt nhalo, offset, size, totalSize;
		MESIOInt *ids, *position;
		MESIOReal *coordinates;

		MESIONodes(mesio, &nhalo, &offset, &size, &totalSize, &ids, &position, &coordinates);
	}

	{ // GET ELEMENTS
		MESIOInt offset, size, totalSize;
		MESIOInt *type, *enodesDist, *enodesData;

		MESIOElements(mesio, &offset, &size, &totalSize, &type, &enodesDist, &enodesData);
	}
	MESIOFinalize();
	MPI_Finalize();
	return 0;
}
```

# License

See the LICENSE file at the root directory.

# Acknowledgment

This work was supported by
 - The Ministry of Education, Youth and Sports from the Large Infrastructures for Research, Experimental Development and Innovations project "e-Infrastructure CZ -- LM2018140",
 - The Ministry of Education, Youth and Sports from the National Programme of Sustainability (NPS II) project "IT4Innovations excellence in science -- LQ1602",
 - The IT4Innovations infrastructure which is supported from the Large Infrastructures for Research, Experimental Development and Innovations project "IT4Innovations National Supercomputing Center -- LM2015070".

