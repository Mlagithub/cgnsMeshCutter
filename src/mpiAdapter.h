#ifndef MPIADAPTER_H
#define MPIADAPTER_H

#ifdef MeshCut_PARALLEL
#endif

#include <mpi.h>
#include <type_traits>

namespace MeshCut
{

typedef MPI_Status MeshCut_Status;
typedef MPI_Aint MeshCut_Aint;
typedef MPI_Request MeshCut_Request;

typedef MPI_Comm MeshCut_Comm;
constexpr MeshCut_Comm MeshCut_COMM_WORLD = MPI_COMM_WORLD;

typedef MPI_Datatype MeshCut_Datatype;
constexpr MeshCut_Datatype MeshCut_DOUBLE = MPI_DOUBLE;
constexpr MeshCut_Datatype MeshCut_INT = MPI_INT;

template<typename DOUBLE, typename std::enable_if<std::is_same<DOUBLE, double>::value, int>::type = 0>
constexpr MPI_Datatype GetMPIDataType()
{
    return MPI_DOUBLE;
}

template <typename FLOAT, typename std::enable_if<std::is_same<FLOAT, float>::value, int>::type = 0>
constexpr MPI_Datatype GetMPIDataType()
{
    return MPI_FLOAT;
}

template<typename INTEGER, typename std::enable_if<std::is_integral<INTEGER>::value, int>::type = 0>
constexpr MPI_Datatype GetMPIDataType()
{
    return std::is_signed<INTEGER>::value ? 
            (
                sizeof(INTEGER) == 1 ?
                MPI_INT8_T
                :(
                sizeof(INTEGER) == 2 ?
                MPI_INT16_T
                :(
                sizeof(INTEGER) == 4 ?
                MPI_INT32_T
                :(
                sizeof(INTEGER) == 8 ?
                MPI_INT64_T
                :
                MPI_DATATYPE_NULL
                )))
            ):(
                sizeof(INTEGER) == 1 ?
                MPI_UINT8_T
                :(
                sizeof(INTEGER) == 2 ?
                MPI_UINT16_T
                :(
                sizeof(INTEGER) == 4 ?
                MPI_UINT32_T
                :(
                sizeof(INTEGER) == 8 ?
                MPI_UINT64_T
                :
                MPI_DATATYPE_NULL
                )))
            );
}

class Mesh;

template<typename Type>
class VolField;

class MPIAdapter
{
private:
    static int _size;
    static int _rank;

public:
    static void Initialize(int *argc, char ***argv);
    static void Finalize();

    static inline int size() { return _size; }
    static inline int rank() { return _rank; }

    static inline bool isMaster() { return _rank == 0; }
    static inline bool isParallel() { return _size > 1; }

    static inline void Barrier() { MPI_Barrier(MPI_COMM_WORLD); }

    template <typename DATATYPE>
    static void Send(const DATATYPE *data, int count, int dest, int tag)
    {
        MPI_Send(data, count, GetMPIDataType<DATATYPE>(), dest, tag, MPI_COMM_WORLD);
    }

    template <typename DATATYPE>
    static void Recv(DATATYPE *data, int count, int src, int tag, MeshCut_Status *status)
    {
        MPI_Recv(data, count, GetMPIDataType<DATATYPE>(), src, tag, MPI_COMM_WORLD, status);
    }

    template <typename DATATYPE>
    static void ISend(const DATATYPE *buf, int count, int dest, int tag, int *request)
    {
        MPI_Isend(buf, count, GetMPIDataType<DATATYPE>(), dest, tag, MPI_COMM_WORLD, request);
    }

    inline static void ISend(const void* buf, int count, MeshCut_Datatype type, int dest, int tag, MeshCut_Comm comm, MeshCut_Request *request)
    {
        MPI_Isend(buf, count, type, dest, tag, comm, request);
    }

    template <typename DATATYPE>
    static void IRecv(DATATYPE *buf, int count, int src, int tag, int *request)
    {
        MPI_Irecv(buf, count, GetMPIDataType<DATATYPE>(), src, tag, MPI_COMM_WORLD, request);
    }

    inline static void IRecv(void* buf, int count, MeshCut_Datatype type, int src, int tag, MeshCut_Comm comm, MeshCut_Request* request)
    {
        MPI_Irecv(buf, count, type, src, tag, comm, request);
    }

    static void WaitAll(int count, int *requests, MeshCut_Status *status);

    static void RequestFree(int *request);

    template <typename DATATYPE>
    static void AllReduceSum(DATATYPE *src, DATATYPE *dst, int count)
    {
        MPI_Allreduce(src, dst, count, GetMPIDataType<DATATYPE>(), MPI_SUM, MPI_COMM_WORLD);
    }

    static void TypeContiguous(const int count, MeshCut_Datatype oldtype, MeshCut_Datatype *newType);
    static void TypeStruct(int count, int *blk_len, MeshCut_Aint *disp, MeshCut_Datatype *types, MeshCut_Datatype *newType);
    static void TypeCommit(MeshCut_Datatype *newType);
    static void TypeFree(MeshCut_Datatype *newType);
};
} // namespace MeshCut

#endif