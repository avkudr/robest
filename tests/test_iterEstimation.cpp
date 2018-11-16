#include "gtest/gtest.h"

#include <random>

#include "LineFitting/LineFitting.hpp"


TEST(iterEstimation, zeroOutliers)
{
    // Estimator initialisation
    robest::RANSAC * ransacSolver = new robest::RANSAC();

    int dataSize = 500;
    double alpha = 0.99;    // 99% success probability
    double gamma = 1;       // 100% of inliers

    int estimNbIter = ransacSolver->estimateNbIter(dataSize, alpha, gamma);
    ASSERT_TRUE(estimNbIter == 1);
}

TEST(iterEstimation, maxAlpha)
{
    // Estimator initialisation
    robest::RANSAC * ransacSolver = new robest::RANSAC();

    int dataSize = 500;
    double alpha = 1;          // 100% success probability
    double gamma = 0.95;       // 95% of inliers

    int estimNbIter = ransacSolver->estimateNbIter(dataSize, alpha, gamma);

    ASSERT_TRUE(estimNbIter == 20000);
}

TEST(iterEstimation, minAlpha)
{
    // Estimator initialisation
    robest::RANSAC * ransacSolver = new robest::RANSAC();

    int dataSize = 500;
    double alpha = -5;         // success probability is incorrect
    double gamma = 0.95;       // 95% of inliers

    int estimNbIter = ransacSolver->estimateNbIter(dataSize, alpha, gamma);

    ASSERT_TRUE(estimNbIter == 20000);
}

TEST(iterEstimation, maxGamma)
{
    // Estimator initialisation
    robest::RANSAC * ransacSolver = new robest::RANSAC();

    int dataSize = 500;
    double alpha = 0.99;         
    double gamma = 5;            // incorrected inlier's ratio

    int estimNbIter = ransacSolver->estimateNbIter(dataSize, alpha, gamma);

    ASSERT_TRUE(estimNbIter == 1);
}

TEST(iterEstimation, zeroInliers)
{
    // Estimator initialisation
    robest::RANSAC * ransacSolver = new robest::RANSAC();

    int dataSize = 500;
    double alpha = 0.99;         
    double gamma = 0;            // incorrected inlier's ratio

    int estimNbIter = ransacSolver->estimateNbIter(dataSize, alpha, gamma);

    ASSERT_TRUE(estimNbIter == 20000);
}

TEST(iterEstimation, zeroIterationInit)
{
    int nbIter = 0;

    std::vector<double> x = {1, 2, 3, 4, 5};
    std::vector<double> y = {1, 2, 3, 4, 5};

    // Define estimation problem
    LineFittingProblem * lineFitting = new LineFittingProblem();
    lineFitting->setData(x,y);

    // Solve
    robest::RANSAC * ransacSolver = new robest::RANSAC();
    ransacSolver->solve(lineFitting, 0.1, nbIter);
}