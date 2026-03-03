#include "mpiAdapter.h"

int MPIAdapter::_size = 1;
int MPIAdapter::_rank = 0;

void MPIAdapter::Initialize(int *argc, char ***argv)
{
    int flag;
    MPI_Initialized(&flag);
    if (!flag)
    {
        MPI_Init(argc, argv);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
}

void MPIAdapter::Finalize()
{
    int flag;
    MPI_Finalized(&flag);
}

void MPIAdapter::RequestFree(MPI_Request *request)
{
    MPI_Request_free(request);
}

void MPIAdapter::TypeContiguous(const int count, MPI_Datatype oldtype, MPI_Datatype *newType)
{
    MPI_Type_contiguous(count, oldtype, newType);
}

void MPIAdapter::TypeStruct(int count, int *blk_len, MPI_Aint *disp, MPI_Datatype *types, MPI_Datatype *newType)
{
    MPI_Type_create_struct(count, blk_len, disp, types, newType);
}

void MPIAdapter::TypeCommit(MPI_Datatype *newType)
{
    MPI_Type_commit(newType);
}

void MPIAdapter::TypeFree(MPI_Datatype *newType)
{
    MPI_Type_free(newType);
}
