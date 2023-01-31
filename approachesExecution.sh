echo approachesExecution started

rm -rf buildTwoDimGreedy
rm -rf buildTwoDimRobin
rm -rf buildOneDim
rm -rf buildThreeDimGreedy
rm -rf buildThreeDimRobin
rm -rf buildTwoDimColor

mkdir buildTwoDimGreedy
mkdir buildTwoDimRobin
mkdir buildOneDim
mkdir buildThreeDimGreedy
mkdir buildThreeDimRobin
mkdir buildTwoDimColor

mkdir approachesOutput

touch approachesOutput/twoDimGreedy.txt
touch approachesOutput/twoDimGreedyOptimum.txt
touch approachesOutput/twoDimRobin.txt
touch approachesOutput/oneDim.txt
touch approachesOutput/threeDimGreedy.txt
touch approachesOutput/threeDimRobin.txt
touch approachesOutput/twoDimColor.txt

cmake . -B buildTwoDimGreedy
cmake . -B buildTwoDimRobin -Dround_robin_distr=1
cmake . -B buildOneDim -Done_dim_tasks=1
cmake . -B buildThreeDimGreedy -Dthree_dim_tasks=1
cmake . -B buildThreeDimRobin -Dround_robin_distr=1
cmake . -B buildTwoDimColor -Dtask_oriented_2d=1

make -C buildTwoDimGreedy -j 20
make -C buildTwoDimRobin -j 20
make -C buildOneDim -j 20
make -C buildThreeDimGreedy -j 20
make -C buildThreeDimRobin -j 20
make -C buildTwoDimColor -j 20

echo twoDimGreedy
echo numThreads: 1 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=1 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt
echo numThreads: 2 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=2 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt
echo numThreads: 4 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=4 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt
echo numThreads: 8 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=8 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt
echo numThreads: 16 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=16 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt

echo twoDimGreedyOptimum
echo numThreads: 4 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=4 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 5 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=5 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 6 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=6 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 7 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=7 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 8 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=8 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 9 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=9 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 10 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=10 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 11 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=11 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 12 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=12 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 13 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=13 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 14 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=14 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 15 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=15 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt
echo numThreads: 16 >> approachesOutput/twoDimGreedyOptimum.txt
OMP_NUM_THREADS=16 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedyOptimum.txt

echo twoDimRobin
echo numThreads: 1 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=1 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt
echo numThreads: 2 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=2 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt
echo numThreads: 4 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=4 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt
echo numThreads: 8 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=8 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt
echo numThreads: 16 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=16 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt

echo oneDim
echo numThreads: 1 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=1 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt
echo numThreads: 2 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=2 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt
echo numThreads: 4 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=4 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt
echo numThreads: 8 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=8 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt
echo numThreads: 16 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=16 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt

echo threeDimGreedy
echo numThreads: 1 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=1 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt
echo numThreads: 2 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=2 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt
echo numThreads: 4 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=4 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt
echo numThreads: 8 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=8 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt
echo numThreads: 16 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=16 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt

echo threeDimRobin
OMP_NUM_THREADS=1 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt
echo numThreads: 2 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=2 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt
echo numThreads: 4 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=4 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt
echo numThreads: 8 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=8 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt
echo numThreads: 16 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=16 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt

echo twoDimColor
echo numThreads: 1 >> approachesOutput/twoDimColor.txt
OMP_NUM_THREADS=1 buildTwoDimColor/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimColor.txt
echo numThreads: 2 >> approachesOutput/twoDimColor.txt
OMP_NUM_THREADS=2 buildTwoDimColor/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimColor.txt
echo numThreads: 4 >> approachesOutput/twoDimColor.txt
OMP_NUM_THREADS=4 buildTwoDimColor/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimColor.txt
echo numThreads: 8 >> approachesOutput/twoDimColor.txt
OMP_NUM_THREADS=8 buildTwoDimColor/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimColor.txt
echo numThreads: 16 >> approachesOutput/twoDimColor.txt
OMP_NUM_THREADS=16 buildTwoDimColor/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimColor.txt