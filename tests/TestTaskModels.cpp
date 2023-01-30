#include <gtest/gtest.h>

#include "data/ParticleContainer.h"
#include "data/Particle.h"
#include "data/ParticleGenerator.h"
#include "omp.h"

#include <vector>

void performTaskModelTest(ParticleContainer& pc){

    std::vector<std::vector<std::vector<std::pair<unsigned long, unsigned long>>>> taskGroup = pc.generate3DTaskModel();

    ASSERT_EQ(taskGroup.size(), 26);

    size_t sumCellInteractions{0};
    for(const auto& task: taskGroup){
        std::unordered_set<unsigned int> cellsVisited{};
        for(const auto& oneThreadPackage: task){

            for(auto& [i,j]: oneThreadPackage){
                ASSERT_FALSE(cellsVisited.contains(i)) << "The cell with index " << i << " got accessed multiple times in the same task";
                ASSERT_FALSE(cellsVisited.contains(j)) << "The cell with index " << j << " got accessed multiple times in the same task";
                cellsVisited.emplace(i);
                cellsVisited.emplace(j);
                sumCellInteractions++;
            }
        }
    }

    size_t cellInteractionReference{0};
    pc.forAllDistinctCellNeighbours([&cellInteractionReference](std::vector<double> &force,
                                                                std::vector<double> &oldForce,
                                                                std::vector<double> &x,
                                                                std::vector<double> &v,
                                                                std::vector<double> &m,
                                                                std::vector<int> &type,
                                                                unsigned long count,
                                                                std::vector<unsigned long> &cell0Items,
                                                                std::vector<unsigned long> &cell1Items,
                                                                std::vector<double> &eps,
                                                                std::vector<double> &sig){
        cellInteractionReference++;
    });

    ASSERT_EQ(cellInteractionReference, sumCellInteractions) << "forAllDistinctCellNeighbours had" << cellInteractionReference << " Neighbouring Cell interactions whereas the taskModel produced " << sumCellInteractions;

}

void performAlternativeTaskModelTest(ParticleContainer& pc){
    std::map<unsigned long, std::unordered_set<unsigned long>> edgesTaken{};

    std::vector<std::vector<std::pair<unsigned long, unsigned long>>> taskGroup = pc.generate2DTaskModelSplitIntoThreads();

    ASSERT_EQ(taskGroup.size(), omp_get_max_threads());

    size_t sumCellInteractions{0};
    for(const auto& task: taskGroup){
        for(auto& [i,j]: task){
            sumCellInteractions++;

            if(edgesTaken.contains(i)){
                ASSERT_FALSE(edgesTaken.at(i).contains(j));
                auto jc = j;
                edgesTaken.at(i).emplace(jc);
            }else{
                auto ic = i;
                edgesTaken.emplace(ic, std::unordered_set<unsigned long>{j});
            }
            if(edgesTaken.contains(j)){
                ASSERT_FALSE(edgesTaken.at(j).contains(i));
                auto ic = i;
                edgesTaken.at(j).emplace(ic);
            }else{
                auto jc = j;
                edgesTaken.emplace(jc, std::unordered_set<unsigned long>{i});
            }
        }
    }

    size_t cellInteractionReference{0};
    pc.forAllDistinctCellNeighbours([&cellInteractionReference](std::vector<double> &force,
                                                                std::vector<double> &oldForce,
                                                                std::vector<double> &x,
                                                                std::vector<double> &v,
                                                                std::vector<double> &m,
                                                                std::vector<int> &type,
                                                                unsigned long count,
                                                                std::vector<unsigned long> &cell0Items,
                                                                std::vector<unsigned long> &cell1Items,
                                                                std::vector<double> &eps,
                                                                std::vector<double> &sig){
        cellInteractionReference++;
    });

    ASSERT_EQ(cellInteractionReference, sumCellInteractions) << "forAllDistinctCellNeighbours had" << cellInteractionReference << " Neighbouring Cell interactions whereas the taskModel produced " << sumCellInteractions;

}


void performTaskOriented2DTest(ParticleContainer& pc){
    std::map<unsigned long, std::unordered_set<unsigned long>> edgesTaken{};
    std::unordered_set<unsigned long> cellsTouchedThisTask{};

    std::vector<std::vector<std::pair<unsigned long, unsigned long>>> taskGroup = pc.generate2DTaskModelColoring();

    ASSERT_EQ(taskGroup.size(), 26);

    size_t sumCellInteractions{0};
    for(const auto& tasks: taskGroup){
        cellsTouchedThisTask.clear();
        for(auto& [i,j]: tasks){
            sumCellInteractions++;

            ASSERT_FALSE(cellsTouchedThisTask.contains(i));
            ASSERT_FALSE(cellsTouchedThisTask.contains(j));
            auto iCopy=i; auto jCopy = j;
            cellsTouchedThisTask.emplace(iCopy);
            cellsTouchedThisTask.emplace(jCopy);

            if(edgesTaken.contains(i)){
                ASSERT_FALSE(edgesTaken.at(i).contains(j));
                auto jc = j;
                edgesTaken.at(i).emplace(jc);
            }else{
                auto ic = i;
                edgesTaken.emplace(ic, std::unordered_set<unsigned long>{j});
            }
            if(edgesTaken.contains(j)){
                ASSERT_FALSE(edgesTaken.at(j).contains(i));
                auto ic = i;
                edgesTaken.at(j).emplace(ic);
            }else{
                auto jc = j;
                edgesTaken.emplace(jc, std::unordered_set<unsigned long>{i});
            }
        }
    }

    size_t cellInteractionReference{0};
    pc.forAllDistinctCellNeighbours([&cellInteractionReference](std::vector<double> &force,
                                                                std::vector<double> &oldForce,
                                                                std::vector<double> &x,
                                                                std::vector<double> &v,
                                                                std::vector<double> &m,
                                                                std::vector<int> &type,
                                                                unsigned long count,
                                                                std::vector<unsigned long> &cell0Items,
                                                                std::vector<unsigned long> &cell1Items,
                                                                std::vector<double> &eps,
                                                                std::vector<double> &sig){
        cellInteractionReference++;
    });

    ASSERT_EQ(cellInteractionReference, sumCellInteractions) << "forAllDistinctCellNeighbours had" << cellInteractionReference << " Neighbouring Cell interactions whereas the taskModel produced " << sumCellInteractions;

}

/**
 * Check if
 * -Tasks created cover all the interactions
 * -the same cell doesn't get used twice within the same job (no potential for race conditions)
 *
 * */
TEST(ParticleContainer, initTaskModel) {
    std::list<Particle> buf;

    Body cub;
    cub.shape = cuboid;
    cub.fixpoint = {0.5, 0.5, 0.5};
    cub.dimensions = {3,3, 3};
    cub.distance = 1;
    cub.mass = 1;
    cub.start_velocity = {0., 0., 0.};
    ParticleGenerator::generateCuboid(cub, 0, buf, 3, 1, 1);
    std::vector bufVec(buf.begin(), buf.end());

    double r_cutoff = 1;
    std::array<double, 3> domainSize{3,3,3};
    unsigned int numCells = domainSize[0]*domainSize[1]*domainSize[2]/ (r_cutoff*r_cutoff*r_cutoff);    //just pick an example that is friendly enough that no rounding stuff messes this up
    ParticleContainer pc(bufVec, domainSize, 1, {}, true);

    performTaskModelTest(pc);
    performAlternativeTaskModelTest(pc);
    performTaskOriented2DTest(pc);

    //second larger test
    cub.dimensions = {100, 100, 100};
    cub.distance = 0.5;
    r_cutoff = 10;
    domainSize = {50, 50, 50};
    buf.clear();
    bufVec.clear();
    ParticleGenerator::generateCuboid(cub, 0, buf, 3, 1,1);
    std::vector bufVec2(buf.begin(), buf.end());

    pc = ParticleContainer(bufVec2, domainSize, r_cutoff, {}, true);

    performTaskModelTest(pc);
    performAlternativeTaskModelTest(pc);
    performTaskOriented2DTest(pc);
}

