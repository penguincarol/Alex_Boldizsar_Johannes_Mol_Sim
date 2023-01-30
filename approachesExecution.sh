echo approachesExecution started

echo numThreads: 1
OMP_NUM_THREADS=1 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1
echo numThreads: 2
OMP_NUM_THREADS=2 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1
echo numThreads: 4
OMP_NUM_THREADS=4 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1
echo numThreads: 8
OMP_NUM_THREADS=8 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1
echo numThreads: 16
OMP_NUM_THREADS=16 buildTwoDimGreedy/MolSim ./input/rti3d_shorterEt.xml -bench file -i 1