
#include "mpiinfo.h"
#include "eslog.hpp"
#include "wrappers/mpi/communication.h"

int mesio::info::mpi::rank = 0;
int mesio::info::mpi::size = 1;
MPI_Comm mesio::info::mpi::comm = MPI_COMM_WORLD;

int mesio::info::mpi::irank = 0;
int mesio::info::mpi::isize = 1;
MPI_Comm mesio::info::mpi::icomm = MPI_COMM_SELF;

int mesio::info::mpi::grank = 0;
int mesio::info::mpi::gsize = 1;
MPI_Comm mesio::info::mpi::gcomm = MPI_COMM_WORLD;

int mesio::info::mpi::threading = 0;

using namespace mesio::info;

void mpi::init(int *argc, char ***argv)
{
	MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &threading);

	mpi::comm = MPI_COMM_WORLD;
	MPI_Comm_rank(mpi::comm, &mpi::rank);
	MPI_Comm_size(mpi::comm, &mpi::size);

	mpi::grank = mpi::rank;
	mpi::gsize = mpi::size;
	mpi::gcomm = mpi::comm;
}

void mpi::init(MPI_Comm comm)
{
	mpi::comm = comm;
	MPI_Comm_rank(mpi::comm, &mpi::rank);
	MPI_Comm_size(mpi::comm, &mpi::size);

	mpi::grank = mpi::rank;
	mpi::gsize = mpi::size;
	mpi::gcomm = mpi::comm;
}

bool mpi::divide(int meshDuplication)
{
	if (meshDuplication == 1) {
		return true;
	}

	if (mesio::info::mpi::size % meshDuplication != 0) {
		return false;
	}

	int color = mpi::rank / (mpi::size / meshDuplication);

	MPI_Comm_split(mpi::gcomm, color, mpi::grank, &mpi::comm);
	MPI_Comm_rank(mpi::comm, &mpi::rank);
	MPI_Comm_size(mpi::comm, &mpi::size);

	MPI_Comm_split(mpi::gcomm, mpi::rank, mpi::grank, &mpi::icomm);
	MPI_Comm_rank(mpi::icomm, &mpi::irank);
	MPI_Comm_size(mpi::icomm, &mpi::isize);

	MPITools::reinit();
	eslog::reinit();

	return true;
}

void mpi::finish()
{
	if (mpi::isize > 1) {
		MPI_Comm_free(&mpi::comm);
		MPI_Comm_free(&mpi::icomm);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}


