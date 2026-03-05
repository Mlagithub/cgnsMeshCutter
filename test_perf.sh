#!/bin/bash
# 性能测试脚本 - 支持多次运行取平均值

export HDF5_DIR=/home/one/app/hdf5-1.12.0/share/cmake/hdf5
source /home/one/projects/cgnsMeshCutter/env.sh

BUILD_DIR=/home/one/projects/cgnsMeshCutter/build
TEST_DIR=/home/one/projects/cgnsMeshCutter/test
CUTTER=$BUILD_DIR/cutter.exe

echo "=========================================="
echo "性能测试: $(date)"
echo "=========================================="

# 清理旧文件
cd $TEST_DIR
rm -f decomposed_mesh_*.cgns residual.dat

# 编译
cd $BUILD_DIR
make -j4 2>&1 | tail -3

# 测试 ParMETIS 分区性能
echo ""
echo "--- ParMETIS 分区性能测试 ---"
cd $TEST_DIR

# 运行多次取平均值
RUNS=5
declare -a TIMES
declare -a MEMS

for i in $(seq 1 $RUNS); do
    echo -n "运行 $i/$RUNS... "
    rm -f decomposed_mesh_*.cgns

    # 使用 /usr/bin/time 获取详细性能数据
    OUTPUT=$(/usr/bin/time -f "ELAPSED:%e\nMEM:%M" mpirun -n 2 $CUTTER -m pipe.cgns -np 2 2>&1)

    ELAPSED=$(echo "$OUTPUT" | grep "ELAPSED:" | cut -d: -f2)
    MEM=$(echo "$OUTPUT" | grep "MEM:" | cut -d: -f2)

    echo "时间: ${ELAPSED}s, 内存: ${MEM}KB"
    TIMES+=($ELAPSED)
    MEMS+=($MEM)
done

# 计算平均值
SUM_TIME=0
SUM_MEM=0
for i in "${!TIMES[@]}"; do
    SUM_TIME=$(echo "$SUM_TIME + ${TIMES[$i]}" | bc)
    SUM_MEM=$(echo "$SUM_MEM + ${MEMS[$i]}" | bc)
done

AVG_TIME=$(echo "scale=4; $SUM_TIME / $RUNS" | bc)
AVG_MEM=$(echo "scale=0; $SUM_MEM / $RUNS" | bc)

# 计算标准差
SUM_SQ=0
for i in "${!TIMES[@]}"; do
    DIFF=$(echo "${TIMES[$i]} - $AVG_TIME" | bc)
    SQ=$(echo "$DIFF * $DIFF" | bc)
    SUM_SQ=$(echo "$SUM_SQ + $SQ" | bc)
done
STD_DEV=$(echo "scale=4; sqrt($SUM_SQ / $RUNS)" | bc)

echo ""
echo "=========================================="
echo "分区性能结果 ($RUNS 次运行):"
echo "  平均时间: ${AVG_TIME}s (标准差: ${STD_DEV}s)"
echo "  平均内存: ${AVG_MEM}KB"
echo "=========================================="

# 测试 ACT 求解器
echo ""
echo "--- ACT 求解器验证 (运行 20 迭代) ---"
rm -f residual.dat
timeout 30 mpirun -n 2 /home/one/app/act/incomACT.exe -c pipe.yaml 2>&1 | grep -E "residual" | head -20

if [ -f residual.dat ]; then
    echo ""
    echo "残差文件已生成，检查收敛性..."
    # 检查残差是否在减小
    FIRST_RES=$(head -2 residual.dat | tail -1)
    LAST_RES=$(tail -1 residual.dat)
    echo "  初始残差: $FIRST_RES"
    echo "  最终残差: $LAST_RES"
fi

echo ""
echo "测试完成: $(date)"