#pragma once

#include <iostream>
#include <vector>
#include <numeric>

#include <mpi.h>

template <typename DOUBLE, typename std::enable_if<std::is_same<DOUBLE, double>::value, int>::type = 0>
constexpr MPI_Datatype GetMPIDataType()
{
    return MPI_DOUBLE;
}

template <typename FLOAT, typename std::enable_if<std::is_same<FLOAT, float>::value, int>::type = 0>
constexpr MPI_Datatype GetMPIDataType()
{
    return MPI_FLOAT;
}

template <typename INTEGER, typename std::enable_if<std::is_integral<INTEGER>::value, int>::type = 0>
constexpr MPI_Datatype GetMPIDataType()
{
    return std::is_signed<INTEGER>::value ? (sizeof(INTEGER) == 1 ? MPI_INT8_T : (sizeof(INTEGER) == 2 ? MPI_INT16_T : (sizeof(INTEGER) == 4 ? MPI_INT32_T : (sizeof(INTEGER) == 8 ? MPI_INT64_T : MPI_DATATYPE_NULL))))
                                          : (sizeof(INTEGER) == 1 ? MPI_UINT8_T : (sizeof(INTEGER) == 2 ? MPI_UINT16_T : (sizeof(INTEGER) == 4 ? MPI_UINT32_T : (sizeof(INTEGER) == 8 ? MPI_UINT64_T : MPI_DATATYPE_NULL))));
}

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
    static inline bool isLast() { return _rank == _size - 1; }
    static inline bool isParallel() { return _size > 1; }

    static inline void Barrier() { MPI_Barrier(MPI_COMM_WORLD); }

    template <typename DATATYPE>
    static void Bcast(DATATYPE *buf, int count, int root)
    {
        MPI_Bcast(buf, count, GetMPIDataType<DATATYPE>(), root, MPI_COMM_WORLD);
    }

    template <typename DATATYPE>
    static void Send(const DATATYPE *data, int count, int dest, int tag)
    {
        MPI_Send(data, count, GetMPIDataType<DATATYPE>(), dest, tag, MPI_COMM_WORLD);
    }

    template <typename DATATYPE>
    static void Recv(DATATYPE *data, int count, int src, int tag, MPI_Status *status)
    {
        MPI_Recv(data, count, GetMPIDataType<DATATYPE>(), src, tag, MPI_COMM_WORLD, status);
    }

    template <typename DATATYPE>
    static void ISend(const DATATYPE *buf, int count, int dest, int tag, MPI_Request *request)
    {
        MPI_Isend(buf, count, GetMPIDataType<DATATYPE>(), dest, tag, MPI_COMM_WORLD, request);
    }

    inline static void ISend(const void *buf, int count, MPI_Datatype type, int dest, int tag, MPI_Comm comm, MPI_Request *request)
    {
        MPI_Isend(buf, count, type, dest, tag, comm, request);
    }

    template <typename DATATYPE>
    static void IRecv(DATATYPE *buf, int count, int src, int tag, MPI_Request *request)
    {
        MPI_Irecv(buf, count, GetMPIDataType<DATATYPE>(), src, tag, MPI_COMM_WORLD, request);
    }

    inline static void IRecv(void *buf, int count, MPI_Datatype type, int src, int tag, MPI_Comm comm, MPI_Request *request)
    {
        MPI_Irecv(buf, count, type, src, tag, comm, request);
    }

    template <typename DATATYPE>
    static void AllReduce(DATATYPE& src, DATATYPE& dst, int count, MPI_Op op)
    {
        if (isParallel())
        {
            MPI_Allreduce(&src, &dst, count, GetMPIDataType<DATATYPE>(), op, MPI_COMM_WORLD);
        }
        else{
            dst = src;
        }
    }
    static void WaitAll(int count, MPI_Request *requests, MPI_Status *status)
    {
        MPI_Waitall(count, requests, status);
    }

    static void RequestFree(MPI_Request *request);

    template <typename DATATYPE>
    static void AllGatherV(DATATYPE* sendBuf, int count, std::vector<DATATYPE>& recvBuf, std::vector<int>& disp)
    {
        int *lenAll = new int[MPIAdapter::size()];
        MPI_Allgather(&count, 1, GetMPIDataType<DATATYPE>(), lenAll, 1, GetMPIDataType<DATATYPE>(), MPI_COMM_WORLD);      
        int lenTotal = std::accumulate(lenAll, lenAll+MPIAdapter::size(), 0);
        recvBuf.assign(lenTotal, 0);
        disp.assign(MPIAdapter::size(),0);
        for(int i=1; i<MPIAdapter::size(); ++i)
        {
            disp[i] = disp[i-1] + lenAll[i-1];
        }
        MPI_Allgatherv(sendBuf, count, GetMPIDataType<DATATYPE>(), recvBuf.data(), lenAll, disp.data(), GetMPIDataType<DATATYPE>(), MPI_COMM_WORLD);

        disp[0] = lenAll[0];
        for(int i=1; i<MPIAdapter::size(); ++i)
        {
            disp[i] = disp[i-1] + lenAll[i]; 
        }
    }

    static void TypeContiguous(const int count, MPI_Datatype oldtype, MPI_Datatype *newType);
    static void TypeStruct(int count, int *blk_len, MPI_Aint *disp, MPI_Datatype *types, MPI_Datatype *newType);
    static void TypeCommit(MPI_Datatype *newType);
    static void TypeFree(MPI_Datatype *newType);
};
