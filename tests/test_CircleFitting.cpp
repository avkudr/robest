#include "gtest/gtest.h"

#include <random>

#include "CircleFitting/CircleFitting.hpp"

//gaussian noise generation
void generateCircleData(
        const double cx,
        const double cy,
        const double r,
        const double noiseVar,
        std::vector<double> & x, std::vector<double> & y)
{
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,noiseVar);

    for(double i = 0 ; i < 3.1415*2.0 ; i += 3.1415/18.0){
        double xnoise = distribution(generator);
        double ynoise = distribution(generator);
        x.push_back(r*cos(double(i)) + cx); //cx
        y.push_back(r*sin(double(i)) + cy); //cy

        if (noiseVar != 0){
            x[i] += xnoise;
            y[i] += ynoise;
        }
    }
}

TEST(CircleFitting, idealCase)
{
    std::vector<double> x;
    std::vector<double> y;

    double cx = 0; // circle center C:(cx,cy)
    double cy = 0;
    double radius = 1; //circle radius
    double noiseVar = 0.0;

    generateCircleData(cx,cy,radius,noiseVar,x,y);

    CircleFittingProblem * circleFitting = new CircleFittingProblem();
    circleFitting->setData(x,y);

    robest::LMedS * solver = new robest::LMedS();
    solver->solve(circleFitting);

    double res_cx,res_cy,res_r;
    circleFitting->getResult(res_cx,res_cy,res_r);

    ASSERT_NEAR(    cx, res_cx, 1.0e-11);
    ASSERT_NEAR(    cy, res_cy, 1.0e-11);
    ASSERT_NEAR(radius,  res_r, 1.0e-11);
}   

TEST(CircleFitting, idealCase2)
{
    std::vector<double> x;
    std::vector<double> y;

    double cx = 24.8; // circle center C:(cx,cy)
    double cy = 8.10;
    double radius = 26.03; //circle radius
    double noiseVar = 0.0;

    generateCircleData(cx,cy,radius,noiseVar,x,y);

    CircleFittingProblem * circleFitting = new CircleFittingProblem();
    circleFitting->setData(x,y);

    double thres = 0.001;
    int nbIter = 20;
    robest::LMedS * solver = new robest::LMedS();
    solver->solve(circleFitting, thres, nbIter);

    double res_cx,res_cy,res_r;
    circleFitting->getResult(res_cx,res_cy,res_r);

    ASSERT_NEAR(    cx, res_cx, 1.0e-11);
    ASSERT_NEAR(    cy, res_cy, 1.0e-11);
    ASSERT_NEAR(radius,  res_r, 1.0e-11);
}

TEST(CircleFitting, smallNoise)
{
    std::vector<double> x;
    std::vector<double> y;

    double cx = 3.552356; // circle center C:(cx,cy)
    double cy = 1.58452;
    double radius = 13.2548; //circle radius
    double noiseVar = 0.001;

    generateCircleData(cx,cy,radius,noiseVar,x,y);

    CircleFittingProblem * circleFitting = new CircleFittingProblem();
    circleFitting->setData(x,y);

    robest::LMedS * solver = new robest::LMedS();
    solver->solve(circleFitting);

    double res_cx,res_cy,res_r;
    circleFitting->getResult(res_cx,res_cy,res_r);

    ASSERT_NEAR(    cx, res_cx, 1.0e-3);
    ASSERT_NEAR(    cy, res_cy, 1.0e-3);
    ASSERT_NEAR(radius,  res_r, 1.0e-3);
}

TEST(CircleFitting, outliers)
{
    std::vector<double> x = {1,0,-1, 0, sqrt(2)/2.0, 24,  8, 26};
    std::vector<double> y = {0,1, 0,-1, sqrt(2)/2.0,  8, 10,  3};

    CircleFittingProblem * circleFitting = new CircleFittingProblem();
    circleFitting->setData(x,y);

    robest::LMedS * solver = new robest::LMedS();
    auto nbIter = solver->calculateIterationsNb(x.size(),0.99,0.45);
    solver->solve(circleFitting, 0.1, nbIter);

    double res_cx,res_cy,res_r;
    circleFitting->getResult(res_cx,res_cy,res_r);

    ASSERT_NEAR( 0.0, res_cx, 1.0e-11);
    ASSERT_NEAR( 0.0, res_cy, 1.0e-11);
    ASSERT_NEAR( 1.0,  res_r, 1.0e-11);
}

TEST(CircleFitting, isDegenerate)
{
    // Generate data
    // y = k*x + b
    std::vector<double> x1 = {0,1,2};
    std::vector<double> y1 = {0,1,2};

    // Define estimation problem
    CircleFittingProblem * circleFitting = new CircleFittingProblem();
    circleFitting->setData(x1,y1);
    ASSERT_TRUE(circleFitting->isDegenerate({0,1,2}));

    std::vector<double> x2 = {0,1,2};
    std::vector<double> y2 = {0,1.001,2};
    circleFitting->setData(x2,y2);
    ASSERT_TRUE(circleFitting->isDegenerate({0,1,2}));

    std::vector<double> x3 = {0,1,5};
    std::vector<double> y3 = {0,1,2};
    circleFitting->setData(x3,y3);
    ASSERT_TRUE(!circleFitting->isDegenerate({0,1,2}));
}
