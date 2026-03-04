#include "mpiAdapter.h"

std::unique_ptr<compi::Environment> MPIAdapter::_env = nullptr;
compi::Context* MPIAdapter::_ctx = nullptr;
int MPIAdapter::_size = 1;
int MPIAdapter::_rank = 0;

void MPIAdapter::Initialize(int *argc, char ***argv)
{
    // Initialize MPI via compi::Environment
    _env = std::make_unique<compi::Environment>(compi::ThreadLevel::Multiple);
    _ctx = &compi::Context::for_comm(MPI_COMM_WORLD);
    
    MPI_Comm_size(MPI_COMM_WORLD, &_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
}

void MPIAdapter::Finalize()
{
    _ctx = nullptr;
    _env.reset();
    _size = 1;
    _rank = 0;
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
