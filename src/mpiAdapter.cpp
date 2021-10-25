#include "mpiAdapter.h"
using namespace std;

namespace MeshCut
{
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
    if (!flag)
    {
        MPI_Finalize();
    }
}

void MPIAdapter::WaitAll(int count, int *requests, MeshCut_Status *status)
{
    MPI_Waitall(count, requests, status);
}

void MPIAdapter::RequestFree(int *request)
{
    MPI_Request_free(request);
}

void MPIAdapter::TypeContiguous(const int count, MeshCut_Datatype oldtype, MeshCut_Datatype *newType)
{
    MPI_Type_contiguous(count, oldtype, newType);
}


void MPIAdapter::TypeStruct(int count, int *blk_len, MeshCut_Aint *disp, MeshCut_Datatype *types, MeshCut_Datatype *newType)
{
    MPI_Type_struct(count, blk_len, disp, types, newType);
}

void MPIAdapter::TypeCommit(MeshCut_Datatype *newType)
{
    MPI_Type_commit(newType);
}

void MPIAdapter::TypeFree(MeshCut_Datatype *newType)
{
    MPI_Type_free(newType);
}
} // namespace MeshCut
