echo approachesExecution started

rm -rf buildTwoDimGreedy
rm -rf buildTwoDimRobin
rm -rf buildOneDim
rm -rf buildThreeDimGreedy
rm -rf buildThreeDimRobin

mkdir buildTwoDimGreedy
mkdir buildTwoDimRobin
mkdir buildOneDim
mkdir buildThreeDimGreedy
mkdir buildThreeDimRobin

mkdir approachesOutput

touch approachesOutput/twoDimGreedy.txt
touch approachesOutput/twoDimRobin.txt
touch approachesOutput/oneDim.txt
touch approachesOutput/threeDimGreedy.txt
touch approachesOutput/threeDimRobin.txt

cmake . -B buildTwoDimGreedy
cmake . -B buildTwoDimRobin -Dround_robin_distr=1
cmake . -B buildOneDim -Done_dim_tasks=1
cmake . -B buildThreeDimGreedy -Dthree_dim_tasks=1
cmake . -B buildThreeDimRobin -Dround_robin_distr=1

make -C buildTwoDimGreedy -j 10
make -C buildTwoDimRobin -j 10
make -C buildOneDim -j 10
make -C buildThreeDimGreedy -j 10
make -C buildThreeDimRobin -j 10

echo testRun to see whether it works. Please ignore: >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=1 buildTwoDimGreedy/MolSim ./input/rti3d_scaled_down.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt

echo numThreads: 1 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=1 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt
echo numThreads: 2 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=2 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt
echo numThreads: 4 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=4 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt
echo numThreads: 8 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=8 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt
echo numThreads: 14 >> approachesOutput/twoDimGreedy.txt
OMP_NUM_THREADS=14 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimGreedy.txt

echo numThreads: 1 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=1 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt
echo numThreads: 2 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=2 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt
echo numThreads: 4 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=4 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt
echo numThreads: 8 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=8 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt
echo numThreads: 14 >> approachesOutput/twoDimRobin.txt
OMP_NUM_THREADS=14 buildTwoDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/twoDimRobin.txt

echo numThreads: 1 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=1 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt
echo numThreads: 2 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=2 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt
echo numThreads: 4 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=4 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt
echo numThreads: 8 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=8 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt
echo numThreads: 14 >> approachesOutput/oneDim.txt
OMP_NUM_THREADS=14 buildOneDim/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/oneDim.txt

echo numThreads: 1 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=1 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt
echo numThreads: 2 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=2 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt
echo numThreads: 4 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=4 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt
echo numThreads: 8 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=8 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt
echo numThreads: 14 >> approachesOutput/threeDimGreedy.txt
OMP_NUM_THREADS=14 buildThreeDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimGreedy.txt

echo numThreads: 1 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=1 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt
echo numThreads: 2 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=2 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt
echo numThreads: 4 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=4 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt
echo numThreads: 8 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=8 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt
echo numThreads: 14 >> approachesOutput/threeDimRobin.txt
OMP_NUM_THREADS=14 buildThreeDimRobin/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1 >> approachesOutput/threeDimRobin.txt