#include "gtest/gtest.h"

#include <random>

#include "LineFitting/LineFitting.hpp"

TEST(iterEstimationDeathTest, wrongInliersRatio)
{
    // Estimator initialisation
    robest::RANSAC * solver = new robest::RANSAC();

    int dataSize = 500;
    double alpha = 0.99;    // 99% success probability

    const auto msg = "Accepted value of inlier's ratio [(]gamma[)] is in range [(]0,1[]]";

    double gamma;
    gamma = 1.1;     // 110% of inliers - error
    EXPECT_DEATH(solver->calculateIterationsNb(dataSize, alpha, gamma), msg);
    gamma = -1;     // negative percentage of inliers - error
    EXPECT_DEATH(solver->calculateIterationsNb(dataSize, alpha, gamma), msg);
    gamma = 0;
    EXPECT_DEATH(solver->calculateIterationsNb(dataSize, alpha, gamma), msg);
}

TEST(iterEstimationDeathTest, wrongSuccessProbability)
{
    // Estimator initialisation
    robest::RANSAC * solver = new robest::RANSAC();
    int dataSize = 500;
    double gamma = 0.95;       // 95% of inliers
    const auto msg = "Accepted value of success probability [(]aplha[)] is in range [(]0,1[)]";

    double alpha;
    alpha =  1.0;          // 100% success probability - error - be more realist
    EXPECT_DEATH(solver->calculateIterationsNb(dataSize, alpha, gamma), msg);
    alpha =  0.0;          //   0% success probability - error - be more optimist
    EXPECT_DEATH(solver->calculateIterationsNb(dataSize, alpha, gamma), msg);
    alpha = -1.0;         //   negative success probability - error - oh common! :)
    EXPECT_DEATH(solver->calculateIterationsNb(dataSize, alpha, gamma), msg);
}

TEST(iterEstimationDeathTest, wrongInputDataSize)
{
    // Estimator initialisation
    robest::RANSAC * solver = new robest::RANSAC();
    
    int dataSize;
    double alpha = 0.99;         
    double gamma = 0.5;   
    const auto msg = "The size of input data must be > 0";

    dataSize = 0;        // 0 points
    EXPECT_DEATH(solver->calculateIterationsNb(dataSize, alpha, gamma), msg);
    dataSize = -15;      // negative amount of points
    EXPECT_DEATH(solver->calculateIterationsNb(dataSize, alpha, gamma), msg);
}

TEST(checkIterNbResults, valuesIterationsNb)
{
    robest::MSAC * solver = new robest::MSAC();
    
    int dataSize = 8;
    double gamma = 0.50;
    double alpha = 0.99;
    EXPECT_EQ(solver->calculateIterationsNb(dataSize, alpha, gamma),1177);

    dataSize = 6;
    gamma    = 0.95;
    alpha    = 0.99;
    EXPECT_EQ(solver->calculateIterationsNb(dataSize, alpha, gamma),4);

    dataSize = 4;
    gamma    = 0.60;
    alpha    = 0.99;
    EXPECT_EQ(solver->calculateIterationsNb(dataSize, alpha, gamma),34);
    
}
