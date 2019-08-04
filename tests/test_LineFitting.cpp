#include "gtest/gtest.h"

#include <random>

#include "LineFitting/LineFitting.hpp"

void generateLineData(
    const double k,
    const double b,
    const double noiseVar,
    const double outliersRatio,
    std::vector<double> &x, std::vector<double> &y)
{
    bool addNoise = noiseVar != 0;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0, 1);
    if (addNoise)
    {
        distribution = std::normal_distribution<double>(0, noiseVar);
    }

    int pointNumber = 50;
    int outlierNumber = (int)(pointNumber*outliersRatio);

    for (int i = 0; i < pointNumber; i++)
    {
        double xnoise = addNoise ? distribution(generator) : 0.0;
        double ynoise = addNoise ? distribution(generator) : 0.0;
        x.push_back((double)i);
        y.push_back(k * x[i] + b);

        if(i < outlierNumber)
        {
            x[i] += xnoise;
            y[i] += ynoise;
        }
    }
}

TEST(LineFitting, idealCase)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x;
    std::vector<double> y;
    double k = 1.0;
    double b = 2.5;
    double noise = 0;
    double outliersRatio = 0.0;
    generateLineData(k, b, noise, outliersRatio, x, y);

    // Define estimation problem
    auto lineFitting = std::make_shared<LineFittingProblem>();
    lineFitting->setData(x, y);

    // Solve
    robest::RANSAC solver;
    solver.solve(lineFitting);

    // Get result
    double res_k, res_b;
    lineFitting->getResult(res_k, res_b);

    ASSERT_NEAR(k, res_k, 1.0e-11);
    ASSERT_NEAR(b, res_b, 1.0e-11);
}

TEST(LineFitting, idealCase2)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x;
    std::vector<double> y;
    double k = 0.7658;
    double b = 4.4785;
    double noise = 0;
    double outliersRatio = 0.0;
    generateLineData(k, b, noise, outliersRatio, x, y);

    // Define estimation problem
    auto lineFitting = std::make_shared<LineFittingProblem>();
    lineFitting->setData(x, y);

    // Solve
    robest::RANSAC solver;
    solver.solve(lineFitting);

    // Get result
    double res_k, res_b;
    lineFitting->getResult(res_k, res_b);

    ASSERT_NEAR(k, res_k, 1.0e-11);
    ASSERT_NEAR(b, res_b, 1.0e-11);
}

TEST(LineFitting, smallNoise)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x;
    std::vector<double> y;
    double k = 1.0;
    double b = 2.5;
    double noise = 0.0001;
    double outliersRatio = 0.15;
    generateLineData(k, b, noise, outliersRatio, x, y);

    // Define estimation problem
    auto lineFitting = std::make_shared<LineFittingProblem>();
    lineFitting->setData(x, y);

    // Solve
    robest::MSAC solver;
    auto nbIter = solver.calculateIterationsNb(lineFitting->getNbMinSamples(),0.99,(1. - outliersRatio));
    solver.solve(lineFitting, 0.1, nbIter);

    // Get result
    double res_k, res_b;
    lineFitting->getResult(res_k, res_b);

    ASSERT_NEAR(k, res_k, 1.0e-3);
    ASSERT_NEAR(b, res_b, 1.0e-3);
}

TEST(LineFitting, outliers)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 5, 14};
    std::vector<double> y = {0, 1, 2, 3, 4, 5, 6, 7, 8, 26, -8};
    double k = 1.0;
    double b = 0.0;

    // Define estimation problem
    auto lineFitting = std::make_shared<LineFittingProblem>();
    lineFitting->setData(x, y);

    // Solve
    double thres = 0.001;
    int nbIter = 20;
    robest::RANSAC solver;
    solver.solve(lineFitting, thres, nbIter);

    // Get result
    double res_k, res_b;
    lineFitting->getResult(res_k, res_b);

    ASSERT_NEAR(k, res_k, 1.0e-11);
    ASSERT_NEAR(b, res_b, 1.0e-11);
}

TEST(LineFitting, almostHalfOutlier)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x = {0, 1, 2, 3, 4, 5, 14, -2, 4, 5};
    std::vector<double> y = {0, 1, 2, 3, 4, 26, -8, 12, 24, 73};
    double k = 1.0;
    double b = 0.0;

    // Define estimation problem
    auto lineFitting = std::make_shared<LineFittingProblem>();
    lineFitting->setData(x, y);

    // Solve
    robest::RANSAC solver;
    int nbIter = solver.calculateIterationsNb(x.size(), 0.99, 0.5);
    solver.solve(lineFitting, 0.1, nbIter);

    // Get result
    double res_k, res_b;
    lineFitting->getResult(res_k, res_b);

    ASSERT_NEAR(k, res_k, 1.0e-11);
    ASSERT_NEAR(b, res_b, 1.0e-11);
}

TEST(LineFitting, noiseCaseMSAC)
{
    // Generate data
    // y = k*x + b
    //std::vector<double> x = {0, 1, 2, 3, 4, 5, 14, -2, 4, 5};
    //std::vector<double> y = {0, 1, 2, 3, 4, 26, -8, 12, 24, 73};

    std::vector<double> x;
    std::vector<double> y;
    double k = 1.0;
    double b = 0.0;
    double noise = 0.0001;
    double outliersRatio = 0.2;
    generateLineData(k, b, noise, outliersRatio, x, y);

    // Define estimation problem
    auto lineFitting = std::make_shared<LineFittingProblem>();
    lineFitting->setData(x, y);

    // Solve
    robest::MSAC solver;
    solver.solve(lineFitting);

    // Get result
    double res_k, res_b;
    lineFitting->getResult(res_k, res_b);

    ASSERT_NEAR(k, res_k, 1.0e-11);
    ASSERT_NEAR(b, res_b, 1.0e-11);
}
