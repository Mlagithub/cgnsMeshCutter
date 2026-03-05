#!/bin/bash

# 并行剖分算法性能测试脚本

BUILD_DIR="../build"
CUTTER="$BUILD_DIR/cutter.exe"

# 测试网格文件
MESH_SMALL="pipe.cgns"           # ~1MB, ~15k cells
MESH_MEDIUM="pipe.256w.cgns"     # ~160MB, ~2.56M cells
MESH_LARGE="pipe.512w.cgns"      # ~320MB, ~5.12M cells

# 结果文件
RESULT_FILE="performance_results.txt"

echo "========================================" | tee $RESULT_FILE
echo "MetisCutter 并行剖分算法性能测试" | tee -a $RESULT_FILE
echo "测试时间: $(date)" | tee -a $RESULT_FILE
echo "========================================" | tee -a $RESULT_FILE
echo "" | tee -a $RESULT_FILE

# 获取网格信息函数
get_mesh_info() {
    local mesh=$1
    local size=$(ls -lh $mesh | awk '{print $5}')
    echo "网格文件: $mesh (大小: $size)"
}

# 运行测试函数
run_test() {
    local mesh=$1
    local np=$2          # 分区数
    local nproc=$3       # MPI进程数
    local test_name=$4   # 测试名称

    echo "----------------------------------------" | tee -a $RESULT_FILE
    echo "测试: $test_name" | tee -a $RESULT_FILE
    echo "网格: $mesh" | tee -a $RESULT_FILE
    echo "分区数: $np, MPI进程数: $nproc" | tee -a $RESULT_FILE

    # 清理旧文件
    rm -f decomposed_mesh_*.cgns

    # 运行测试并计时
    local start_time=$(date +%s.%N)

    if [ $nproc -eq 1 ]; then
        # 串行测试
        $CUTTER -m $mesh -np $np 2>&1 | tee -a $RESULT_FILE
    else
        # 并行测试
        mpirun -n $nproc $CUTTER -m $mesh -np $np 2>&1 | tee -a $RESULT_FILE
    fi

    local end_time=$(date +%s.%N)
    local elapsed=$(echo "$end_time - $start_time" | bc)

    echo "执行时间: ${elapsed} 秒" | tee -a $RESULT_FILE

    # 统计生成的文件
    local file_count=$(ls decomposed_mesh_*.cgns 2>/dev/null | wc -l)
    echo "生成文件数: $file_count" | tee -a $RESULT_FILE
    echo "" | tee -a $RESULT_FILE

    # 清理
    rm -f decomposed_mesh_*.cgns
}

# ========== 测试1: 小网格性能测试 ==========
echo "========== 测试1: 小网格性能测试 ==========" | tee -a $RESULT_FILE
get_mesh_info $MESH_SMALL | tee -a $RESULT_FILE
echo "" | tee -a $RESULT_FILE

run_test $MESH_SMALL 2 1 "串行-2分区"
run_test $MESH_SMALL 2 2 "并行-2分区-2进程"
run_test $MESH_SMALL 4 2 "并行-4分区-2进程"
run_test $MESH_SMALL 4 4 "并行-4分区-4进程"
run_test $MESH_SMALL 8 4 "并行-8分区-4进程"

# ========== 测试2: 中等网格性能测试 ==========
echo "========== 测试2: 中等网格性能测试 ==========" | tee -a $RESULT_FILE
get_mesh_info $MESH_MEDIUM | tee -a $RESULT_FILE
echo "" | tee -a $RESULT_FILE

run_test $MESH_MEDIUM 2 1 "串行-2分区"
run_test $MESH_MEDIUM 2 2 "并行-2分区-2进程"
run_test $MESH_MEDIUM 4 2 "并行-4分区-2进程"
run_test $MESH_MEDIUM 4 4 "并行-4分区-4进程"
run_test $MESH_MEDIUM 8 4 "并行-8分区-4进程"

# ========== 测试3: 大网格性能测试 ==========
echo "========== 测试3: 大网格性能测试 ==========" | tee -a $RESULT_FILE
get_mesh_info $MESH_LARGE | tee -a $RESULT_FILE
echo "" | tee -a $RESULT_FILE

run_test $MESH_LARGE 2 2 "并行-2分区-2进程"
run_test $MESH_LARGE 4 2 "并行-4分区-2进程"
run_test $MESH_LARGE 4 4 "并行-4分区-4进程"
run_test $MESH_LARGE 8 4 "并行-8分区-4进程"

echo "========================================" | tee -a $RESULT_FILE
echo "性能测试完成" | tee -a $RESULT_FILE
echo "结果保存在: $RESULT_FILE" | tee -a $RESULT_FILE
echo "========================================" | tee -a $RESULT_FILE