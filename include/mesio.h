
#ifndef _MESIO_H_
#define _MESIO_H_

#include "mpi.h"

/*-----------------------------------------------------------------------------
 Set data-types used in MESIO

 Possible values for MESIO_INT_WIDTH:
    32: use 32 bit signed integer
    64: use 64 bit signed integer

 Possible values for MESIO_REAL_WIDTH:
    64: MESIO supports only 64 bit real values
------------------------------------------------------------------------------*/
#ifndef MESIO_INT_WIDTH
#define MESIO_INT_WIDTH 32
#endif

#ifndef MESIO_REAL_WIDTH
#define MESIO_REAL_WIDTH 64
#endif

#if MESIO_INT_WIDTH == 32
    typedef int MESIOInt;
#elif MESIO_INT_WIDTH == 64
    typedef long MESIOInt;
#else
#error "Incorrect user-supplied value of MESIO_INT_WIDTH"
#endif

#if MESIO_REAL_WIDTH == 64
    typedef double MESIOReal;
#else
#error "Incorrect user-supplied value of MESIO_REAL_WIDTH"
#endif


typedef enum {
    MESIO_ANSYS,
    MESIO_ENSIGHT,
    MESIO_VTK_LEGACY,
    MESIO_XDMF,
} MESIOFormat;

typedef enum {
    MESIO_NONE,
    MESIO_METIS,
    MESIO_PARMETIS,
    MESIO_PTSCOTCH,
    MESIO_HILBERT_CURVE
} MESIODecomposer;

typedef enum {
    POINT1, // 0

    // without mid-points
    LINE2, // 1
    TRIANGLE3, // 2
    SQUARE4, // 3
    TETRA4, // 4
    PYRAMID5, // 5
    PRISMA6, // 6
    HEXA8, // 7

    // with mid-points
    LINE3, // 8
    TRIANGLE6, // 9
    SQUARE8, // 10
    TETRA10, // 11
    PYRAMID13, // 12
    PRISMA15, // 13
    HEXA20, // 14

    // element with unknown number of nodes
    NOT_SUPPORTED,

    SIZE
} MESIOElementType;

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
 Definitions of internal structures used in MESIO
------------------------------------------------------------------------------*/
typedef struct MESIOData* MESIO;

/*-----------------------------------------------------------------------------
 Functions for manipulating with MESIO internal structures
------------------------------------------------------------------------------*/

/// Initialize the MESIO library
/**
 * MESIO needs to be initialized before calling any other function.
 * At the end of the program MESIOFinalize should be called.
 *
 * @param comm an MPI communication
 * @param verbosity a verbosity level
 */
void MESIOInit(
    MPI_Comm        comm,
    int             verbosity
);

/// Finalize the MESIO library
/**
 * This function destroy all internal parameters.
 */
void MESIOFinalize();


/// Load an input database and build a mesh
/**
 * This function load an input database from the provided 'path'.
 * The 'path' parameter should point to a database in provided 'format'
 * (see wiki for path requirements for each format). After calling this function,
 * a mesh is stored in the 'mesio' handler. Then, mesh information can be queried
 * by other functions from this handler. One can influence distribution of elements
 * among available processes by 'decomposer' parameter. Elements are distributed
 * according to a selected element node coordinates and Hilbert's space filling
 * curve if 'decomposer=MESIO_NONE'. If Mesio 'decomposer=MESIO_HILBERT_CURVE', then
 * coordinates for Hilbert's curve are computed from elements centers.
 * Mesio also provides second level decomposition (decomposition of elements
 * within each MPI process) by parameter 'domains'. The second level decomposition
 *  can be skipped by setting 'domains=0'.
 *
 * @param mesio handler
 * @param format format of an input database
 * @param path path to an input database
 * @param decomposer internally used parallel decomposer
 * @param domains
 */
void MESIOLoad(
    MESIO*          mesio,
    MESIOFormat     format,
    const char*     path,
    MESIODecomposer decomposer,
    int             domains
);

/// Returns nodes of the loaded mesh
/**
 * Returns nodes description for each MPI processes.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio handler with the loaded mesh
 * @param nhalo number of halo nodes - nodes that are also stored to processes with lower MPI ranks
 * @param offset number of unique nodes stored on processes with lower MPI ranks
 * @param size number of unique nodes that are stored on this process (and possibly by processes with higher MPI rank)
 * @param totalSize number of unique nodes of the loaded mesh
 * @param ids array of node IDs with the length 'nhalo+size'
 * @param position array of node position (a global unique offset) with the length 'nhalo+size'
 * @param coordinates array of node coordinates with the length 'nhalo+size'
 */
void MESIONodes(
    MESIO           mesio,
    MESIOInt*       nhalo,
    MESIOInt*       offset,
    MESIOInt*       size,
    MESIOInt*       totalSize,
    MESIOInt**      ids,
    MESIOInt**      position,
    MESIOReal**     coordinates
);

/// Returns nodes ranks of the loaded mesh
/**
 * Returns ranks of each node. The rank denotes MPI rank of process
 * where the node is stored. Nodes are stored on more MPI processes
 * since nodes are duplicated to all processes where they occurs
 * in any element. Ranks are returned by as two arrays. The first array
 * returns distribution of ranks that are stored in 'rankData' array.
 * The ranks of the first node start at 'rankDistribution[0]' and end
 * at 'rankDistribution[1]'. The ranks of the second node start
 * at 'rankDistribution[1]' and end at 'rankDistribution[2]', and so on.
 * There is always a lower rank for all nodes with offset lower that 'nhalo'.
 * For the rest of nodes the first rank is always the rank of current process.
 * Hence, for each 'i<nhalo' 'rankData[rankDistribution[i]]' is lower than
 * the rank of process that call this function. And, for each 'j>=nhalo'
 * 'rankData[rankDistribution[j]]' is equal to the rank of process that
 * call this function.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio handler with the loaded mesh
 * @param rankDistribution array of length 'nhalo+size+1' with the distribution of ranks
 * @param rankData array of length 'rankDistribution[nhalo+size]' with the ranks data
 */

void MESIONodesRanks(
    MESIO           mesio,
    MESIOInt**      rankDistribution,
    int**           rankData
);

/// Returns nodes domains of the loaded mesh
/**
 * Returns domains of each node. The behavior of this function
 * is the same as the behavior of function 'MESIONodesRanks'.
 *
 * @param mesio handler with the loaded mesh
 * @param domainDistribution array of length 'nhalo+size+1' with the distribution of domains
 * @param domainData array of length 'domainDistribution[nhalo+size]' with the domains data
 */
void MESIONodesDomains(
    MESIO           mesio,
    MESIOInt**      domainDistribution,
    MESIOInt**      domainData
);

/// Returns the node to elements map
/**
 * Returns the list of elements where a node occurs. Elements are returned
 * as global offset to an element. The first array returns distribution of elements
 * that are stored in 'elementData' array. The elements of the first
 * node start at 'elementDistribution[0]' and end at 'elementDistribution[1]'.
 * The elements of the second node start at 'elementDistribution[1]'
 * and end at 'elementDistribution[2]', and so on.
 *
 * @param mesio handler with the loaded mesh
 * @param elementDistribution array of length 'size + 1'
 * @param elementData array of length 'elementDistribution[size]' with the node elements
 */
void MESIONodeToElements(
    MESIO           mesio,
    MESIOInt**      elementDistribution,
    MESIOInt**      elementData
);

/// Returns elements of the loaded mesh
/**
 * Returns elements description for each MPI processes. Elements are stored
 * uniquely on each process.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio handler with the loaded mesh
 * @param offset number of elements stored on processes with lower MPI ranks
 * @param size number of element that are stored on this process
 * @param totalSize number of elements of the loaded mesh
 * @param type array of element type with the length 'size'
 * @param enodesDistribution array of element nodes distribution with the length 'size + 1'
 * @param enodesData array of element nodes with the length 'enodesDistribution[size]'
 */
void MESIOElements(
    MESIO           mesio,
    MESIOInt*       offset,
    MESIOInt*       size,
    MESIOInt*       totalSize,
    MESIOInt**      type,
    MESIOInt**      enodesDistribution,
    MESIOInt**      enodesData
);

/// Returns elements domains
/**
 * Returns the domain index for each element.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio handler with the loaded mesh
 * @param domains array with length 'size'
 */
void MESIOElementsDomains(
    MESIO           mesio,
	MESIOInt**      domains
);

/// Returns elements material
/**
 * Returns the material index for each element.
 *
 * @param mesio handler with the loaded mesh
 * @param material array with length 'size'
 */
void MESIOElementsMaterials(
    MESIO           mesio,
    int**           material
);

/// Returns elements body
/**
 * Returns the body index for each element.
 *
 * @param mesio handler with the loaded mesh
 * @param bodies total number of bodies
 * @param bodies array with length 'size'
 */
void MESIOElementsBodies(
    MESIO           mesio,
    MESIOInt*       bodies,
    int**           body
);

/// Returns elements face neighbors
/**
 * Returns face neighbors for each element. Neighbors are returned
 * as global offset to an element that shared a particular face.
 * Faces order can be listed by function 'MESIOElementFaceList'.
 * If '-1' is returned as a neighbor, then the face is
 * a face on surface. The first array returns distribution of faces
 * that are stored in 'neighborData' array. The faces of the first
 * element start at 'neighborDistribution[0]' and end
 * at 'neighborDistribution[1]'. The faces of the second element start
 * at 'neighborDistribution[1]' and end at 'neighborDistribution[2]',
 * and so on.
 *
 * @param mesio handler with the loaded mesh
 * @param neighborDistribution array of length 'size + 1'
 * @param neighborData array of length 'neighborDistribution[size]' with the face neighbors
 */
void MESIOElementsFaceNeighbors(
    MESIO           mesio,
    MESIOInt**      neighborDistribution,
    MESIOInt**      neighborData
);

/// Returns elements edge neighbors
/**
 * Returns edge neighbors for each element. Neighbors are returned
 * as global offset to an element that shared a particular edge.
 * Edges order can be listed by function 'MESIOElementEdgeList'.
 * The first array returns distribution of edges that are stored
 * in 'neighborData' array. The edges of the first element start
 * at 'neighborDistribution[0]' and end at 'neighborDistribution[1]'.
 * The edges of the second element start at 'neighborDistribution[1]'
 * and end at 'neighborDistribution[2]', and so on. Since there
 * can be more neighboring elements for each edge, the first
 * number denotes the number of neighboring elements (e.g., 2 17 73
 * denotes that there are 2 neighboring elements on a particular
 * edge - elements 17 and 73).
 *
 * @param mesio handler with the loaded mesh
 * @param neighborDistribution array of length 'size + 1'
 * @param neighborData array of length 'neighborDistribution[size]' with the edge neighbors
 */
void MESIOElementsEdgeNeighbors(
    MESIO           mesio,
    MESIOInt**      neighborDistribution,
    MESIOInt**      neighborData
);

/// Returns the list of faces
/**
 * Returns the list of nodes that form a face of an element.
 * The first face is composed of nodes stored form
 * faceDistribution[0] to faceDistribution[1].
 *
 * @param mesio handler with the loaded mesh
 * @param element element type
 * @faceDistribution face distribution
 * @faceData node indices of a particular face
 */
void MESIOElementFaceList(
    MESIO               mesio,
    MESIOElementType    element,
    MESIOInt**          faceDistribution,
    MESIOInt**          faceData
);

/// Returns the list of edges
/**
 * Returns the list of nodes that form an edge of an element.
 * The first edge is composed of nodes stored form
 * edgeDistribution[0] to edgeDistribution[1].
 *
 * @param mesio handler with the loaded mesh
 * @param element element type
 * @edgeDistribution edge distribution
 * @edgeData node indices of a particular edge
 */
void MESIOElementEdgeList(
    MESIO               mesio,
    MESIOElementType    element,
    MESIOInt**          edgeDistribution,
    MESIOInt**          edgeData
);

/// Returns total number of element of a given type
/**
 * Returns total number of element of a given type
 * and the number of elements stored on processes
 * with lower rank (offset).
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @param etype element type
 * @param offset number of element of a given type stored on processes with lower rank
 * @param totalSize number of element of a given type in the mesh
 */
void MESIOElementsCounters(
    MESIO           mesio,
    MESIOInt        etype,
    MESIOInt*       offset,
    MESIOInt*       totalSize
);

/// Returns number of elements regions
/**
 * Returns number of elements regions.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @return number of elements regions
 */
int MESIOElementsRegions(
    MESIO           mesio
);

/// Returns description of a region with a given number
/**
 * Returns description of a region with a given number.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @param region region index
 * @param name name of the region
 * @param size number of elements in a particular region
 * @param elements local offset of elements in a particular region
 */
void MESIOElementsRegion(
    MESIO           mesio,
    MESIOInt        region,
    const char**    name,
    MESIOInt*       size,
    MESIOInt**      elements
);

/// Returns list of nodes of a given region
/**
 * Returns list of nodes of a given region. The function
 * has the same behavior as function 'MESIONodes' except that
 * only nodes of a given region are returned.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @param region region index
 * @param nhalo number of halo nodes - nodes that are also stored to processes with lower MPI ranks
 * @param offset number of unique nodes stored on processes with lower MPI ranks
 * @param size number of unique nodes that are stored on this process (and possibly by processes with higher MPI rank)
 * @param totalSize number of unique nodes of a given region of the loaded mesh
 * @param nodes array of node local offsets with the length 'nhalo+size'
 * @param position array of node position (a global unique offset within a region) with the length 'nhalo+size'
 */
void MESIOElementsRegionNodes(
    MESIO           mesio,
    MESIOInt        region,
    MESIOInt*       nhalo,
    MESIOInt*       offset,
    MESIOInt*       size,
    MESIOInt*       totalSize,
    MESIOInt**      nodes,
	MESIOInt**      position
);

/// Returns total number of element of a given type in a given region
/**
 * Returns total number of element of a given type in a given region
 * and the number of elements stored on processes with lower rank (offset).
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @param region region index
 * @param etype element type
 * @param offset number of element of a given type stored on processes with lower rank
 * @param totalSize number of element of a given type in the mesh
 */
void MESIOElementsRegionCounters(
    MESIO           mesio,
    MESIOInt        region,
    MESIOInt        etype,
    MESIOInt*       offset,
    MESIOInt*       totalSize
);

/// Returns number of boundary regions
/**
 * Returns number of boundary regions.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @return number of boundary regions
 */
int MESIOBoundaryRegions(
    MESIO           mesio
);

/// Returns description of a boundary region with a given number
/**
 * Returns description of a boundary region with a given number.
 * Elements of a boundary region are returned as a list of its nodes
 * since boundary elements are not part of the rest elements.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @param region region index
 * @param name name of the region
 * @param dimension dimension of a given element (0 - points, 1 - edges, 2 - faces)
 * @param size number of elements in a particular region
 * @param type array of element types
 * @param parent array of parents elements
 * @param elementDistribution distribution of elements
 * @param elementData array of elements nodes
 */
void MESIOBoundaryRegion(
    MESIO           mesio,
    MESIOInt        region,
    const char**    name,
    MESIOInt*       dimension,
    MESIOInt*       size,
    MESIOInt**      type,
    MESIOInt**      parent,
    MESIOInt**      elementDistribution,
    MESIOInt**      elementData
);

/// Returns list of nodes of a given boundary region
/**
 * Returns list of nodes of a given boundary region. The function
 * has the same behavior as function 'MESIONodes' except that
 * only nodes of a given boundary region are returned.
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @param region region index
 * @param nhalo number of halo nodes - nodes that are also stored to processes with lower MPI ranks
 * @param offset number of unique nodes stored on processes with lower MPI ranks
 * @param size number of unique nodes that are stored on this process (and possibly by processes with higher MPI rank)
 * @param totalSize number of unique nodes of a given region of the loaded mesh
 * @param nodes array of node local offsets with the length 'nhalo+size'
 * @param position array of node position (a global unique offset within a region) with the length 'nhalo+size'
 */
void MESIOBoundaryRegionNodes(
    MESIO           mesio,
    MESIOInt        region,
    MESIOInt*       nhalo,
    MESIOInt*       offset,
    MESIOInt*       size,
    MESIOInt*       totalSize,
    MESIOInt**      nodes,
    MESIOInt**      position
);

/// Returns total number of element of a given type in a given boundary region
/**
 * Returns total number of element of a given type in a given boundary region
 * and the number of elements stored on processes with lower rank (offset).
 *
 * Please run 'test.mesio' to see a possible output of a simple mesh database.
 *
 * @param mesio mesio handler with the loaded mesh
 * @param region region index
 * @param etype element type
 * @param offset number of element of a given type stored on processes with lower rank
 * @param totalSize number of element of a given type in the mesh
 */
void MESIOBoundaryRegionCounters(
    MESIO           mesio,
    MESIOInt        region,
    MESIOInt        etype,
    MESIOInt*       offset,
    MESIOInt*       totalSize
);

/*-----------------------------------------------------------------------------
 Destroy an arbitrary internal structure
------------------------------------------------------------------------------*/

/// Destroy a data holder
/**
 * Destroy all data associated to a data holder
 *
 * @param ptr an address to a data holder
 */
void MESIODestroy(
    void*        ptr
);

#ifdef __cplusplus
}
#endif

#endif /* _MESIO_H_ */
