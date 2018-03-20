#include "gtest/gtest.h"

#include <random>

#include "LineFitting/LineFitting.hpp"

void generateLineData(
        const double k,
        const double b,
        const double noiseVar,
        std::vector<double> & x, std::vector<double> & y)
{
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,noiseVar);

    for (int i = 0; i < 50; i++)
    {
        double xnoise = distribution(generator);
        double ynoise = distribution(generator);
        x.push_back((double)i);
        y.push_back(k*x[i] + b);

        if (noiseVar != 0){
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
    generateLineData(k,b,noise,x,y);

    // Define estimation problem
    LineFittingProblem * lineFitting = new LineFittingProblem();
    lineFitting->setData(x,y);

    // Solve
    robest::RANSAC * ransacSolver = new robest::RANSAC();
    ransacSolver->solve(lineFitting);

    // Get result
    double res_k,res_b;
    lineFitting->getResult(res_k,res_b);

    ASSERT_NEAR(k,res_k, 1.0e-11);
    ASSERT_NEAR(b,res_b, 1.0e-11);
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
    generateLineData(k,b,noise,x,y);

    // Define estimation problem
    LineFittingProblem * lineFitting = new LineFittingProblem();
    lineFitting->setData(x,y);

    // Solve
    robest::RANSAC * ransacSolver = new robest::RANSAC();
    ransacSolver->solve(lineFitting);

    // Get result
    double res_k,res_b;
    lineFitting->getResult(res_k,res_b);

    ASSERT_NEAR(k,res_k, 1.0e-11);
    ASSERT_NEAR(b,res_b, 1.0e-11);
}

TEST(LineFitting, smallNoise)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x;
    std::vector<double> y;
    double k = 1.0;
    double b = 2.5;
    double noise = 0.001;
    generateLineData(k,b,noise,x,y);

    // Define estimation problem
    LineFittingProblem * lineFitting = new LineFittingProblem();
    lineFitting->setData(x,y);

    // Solve
    robest::RANSAC * ransacSolver = new robest::RANSAC();
    ransacSolver->solve(lineFitting);

    // Get result
    double res_k,res_b;
    lineFitting->getResult(res_k,res_b);

    ASSERT_NEAR(k,res_k, 1.0e-3);
    ASSERT_NEAR(b,res_b, 1.0e-3);
}

TEST(LineFitting, outliers)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x = {0,1,2,3,4,5,6,7,8,  5, 14};
    std::vector<double> y = {0,1,2,3,4,5,6,7,8, 26, -8};
    double k = 1.0;
    double b = 0.0;

    // Define estimation problem
    LineFittingProblem * lineFitting = new LineFittingProblem();
    lineFitting->setData(x,y);

    // Solve
    double thres = 0.001;
    int nbIter = 20;
    robest::RANSAC * ransacSolver = new robest::RANSAC();
    ransacSolver->solve(lineFitting, thres, nbIter);

    // Get result
    double res_k,res_b;
    lineFitting->getResult(res_k,res_b);

    ASSERT_NEAR(k,res_k, 1.0e-11);
    ASSERT_NEAR(b,res_b, 1.0e-11);
}

TEST(LineFitting, almostHalfOutlier)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x = {0,1,2,3,4, 5,14,-2, 4, 5};
    std::vector<double> y = {0,1,2,3,4,26,-8,12,24,73};
    double k = 1.0;
    double b = 0.0;

    // Define estimation problem
    LineFittingProblem * lineFitting = new LineFittingProblem();
    lineFitting->setData(x,y);

    // Solve
    robest::RANSAC * ransacSolver = new robest::RANSAC();
    ransacSolver->solve(lineFitting);

    // Get result
    double res_k,res_b;
    lineFitting->getResult(res_k,res_b);

    ASSERT_NEAR(k,res_k, 1.0e-11);
    ASSERT_NEAR(b,res_b, 1.0e-11);
}
