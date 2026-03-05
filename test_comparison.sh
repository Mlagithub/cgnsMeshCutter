#!/bin/bash
# 串并行剖分性能对比测试

export HDF5_DIR=/home/one/app/hdf5-1.12.0/share/cmake/hdf5
source /home/one/projects/cgnsMeshCutter/env.sh

BUILD_DIR=/home/one/projects/cgnsMeshCutter/build
TEST_DIR=/home/one/projects/cgnsMeshCutter/test
CUTTER=$BUILD_DIR/cutter.exe

MESH_FILE="pipe.256w.cgns"
NPARTS=(2 4 8)

echo "=========================================="
echo "串并行剖分性能对比测试"
echo "网格文件: $MESH_FILE"
echo "日期: $(date)"
echo "=========================================="

cd $TEST_DIR

# 清理
rm -f decomposed_mesh_*.cgns

for np in "${NPARTS[@]}"; do
    echo ""
    echo "=========================================="
    echo "分区数: $np"
    echo "=========================================="

    # 串行测试 (METIS)
    echo ""
    echo "--- 串行 METIS 测试 ---"
    rm -f decomposed_mesh_*.cgns

    /usr/bin/time -v $CUTTER -m $MESH_FILE -np $np 2>&1 | tee /tmp/metis_output.log | grep -E "(METIS|write section|Elapsed|Maximum resident)"

    # 并行测试 (ParMETIS)
    echo ""
    echo "--- 并行 ParMETIS 测试 ($np 进程) ---"
    rm -f decomposed_mesh_*.cgns

    /usr/bin/time -v mpirun -n $np $CUTTER -m $MESH_FILE -np $np 2>&1 | tee /tmp/parmetis_output.log | grep -E "(ParMETIS|write section|Elapsed|Maximum resident)"

    echo ""
done

echo ""
echo "=========================================="
echo "测试完成: $(date)"
echo "=========================================="