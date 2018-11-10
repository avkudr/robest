#include "gtest/gtest.h"

#include <random>

#include "SphereFitting/SphereFitting.hpp"

// generate data with gaussian noise (optional)
void generateSphereData(
        const double cx,
        const double cy,
        const double cz,
        const double r,
        const double noiseVar,
        std::vector<double> & x, std::vector<double> & y, std::vector<double> & z)
{
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,noiseVar);

    for(double i = 0 ; i < 3.1415*2.0 ; i += 3.1415/18.0){
        double xnoise = distribution(generator);
        double ynoise = distribution(generator);
        double znoise = distribution(generator);
        x.push_back(r*cos(double(i))*sin(double(i)) + cx); //cx
        y.push_back(r*sin(double(i))*sin(double(i)) + cy); //cy
        z.push_back(r*cos(double(i)) + cz); //cz

        if (noiseVar != 0){
            x[i] += xnoise;
            y[i] += ynoise;
            z[i] += znoise;
        }
    }
}

TEST(SphereFitting, idealCase)
{
    // Define estimation problem
    SphereFittingProblem * sphereFitting = new SphereFittingProblem();

    // Define dataset
    std::vector<double> x = {0,  2,  0,   0};
    std::vector<double> y = {0,  0,  2,   0};
    std::vector<double> z = {2,  0,  0,  -2};

    double cx = 0; // sphere center C:(cx,cy, cz)
    double cy = 0;
    double cz = 0;
    double radius = 2; //sphere radius
    double noiseVar = 0.0;

    sphereFitting->setData(x,y,z);

    // Solver init
    robest::MSAC * MSACsolver = new robest::MSAC();
    MSACsolver->solve(sphereFitting);

    // Get results
    double res_cx,res_cy, res_cz,res_r;
    sphereFitting->getResult(res_cx, res_cy, res_cz, res_r);

    ASSERT_NEAR(    cx, res_cx, 1.0e-11);
    ASSERT_NEAR(    cy, res_cy, 1.0e-11);
    ASSERT_NEAR(    cz, res_cz, 1.0e-11);
    ASSERT_NEAR(radius,  res_r, 1.0e-11);
}   

TEST(SphereFitting, idealCase2)
{
    // Define estimation problem
    SphereFittingProblem * sphereFitting = new SphereFittingProblem();

    // Define dataset
    std::vector<double> x = {   8,   12,     4,     8,     8,     8};
    std::vector<double> y = {  10,    6,     6,     2,     6,     6};
    std::vector<double> z = { -11,  -11,   -11,   -11,   -15,    -7};

    double cx = 8.0; // sphere center C:(cx,cy,cz)
    double cy = 6.0;
    double cz = -11.0;
    double radius = 4.0; //sphere radius
    double noiseVar = 0.0;

    sphereFitting->setData(x,y,z);

    // Solver init
    robest::MSAC * MSACsolver = new robest::MSAC();
    MSACsolver->solve(sphereFitting);

    // Get results
    double res_cx,res_cy, res_cz,res_r;
    sphereFitting->getResult(res_cx, res_cy, res_cz, res_r);
     
    ASSERT_NEAR(    cx, res_cx, 1.0e-11);
    ASSERT_NEAR(    cy, res_cy, 1.0e-11);
    ASSERT_NEAR(    cz, res_cz, 1.0e-11);
    ASSERT_NEAR(radius,  res_r, 1.0e-11);
}   

TEST(SphereFitting, smallNoise)
{
    // Define estimation problem
    SphereFittingProblem * sphereFitting = new SphereFittingProblem();

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    double cx = 5.23542; //sphere center C:(cx,cy,cz)
    double cy = 3.35662;
    double cz = -1.96731;
    double radius = 4.54854; //sphere radius
    double noiseVar = 0.001;

    // generate dataset
    generateSphereData(cx,cy,cz,radius,noiseVar,x,y,z);

    sphereFitting->setData(x,y,z);

    // Solver init
    robest::MSAC * MSACsolver = new robest::MSAC();
    MSACsolver->solve(sphereFitting);

    // Get results
    double res_cx,res_cy, res_cz,res_r;
    sphereFitting->getResult(res_cx, res_cy, res_cz, res_r);
     
    ASSERT_NEAR(    cx, res_cx, 1.0e-3);
    ASSERT_NEAR(    cy, res_cy, 1.0e-3);
    ASSERT_NEAR(    cz, res_cz, 1.0e-3);
    ASSERT_NEAR(radius,  res_r, 1.0e-3);
}

TEST(SphereFitting, outliers)
{
    // Define estimation problem
    SphereFittingProblem * sphereFitting = new SphereFittingProblem();

    // Define dataset
    std::vector<double> x = {3, 8,  -2, 3,  3,  3,  16, -24, 0};
    std::vector<double> y = {3, -2, -2, -7, -2, -2, 12, 16, -15, 0};
    std::vector<double> z = {7, 7,  7,  7,  2,  12, -2, 3,  0};

    double cx = 3; // sphere center C:(cx,cy,cz)
    double cy = -2;
    double cz = 7;
    double radius = 5; //sphere radius;

    sphereFitting->setData(x,y,z);

    // Solver init
    robest::MSAC * MSACsolver = new robest::MSAC();
    MSACsolver->solve(sphereFitting);

    // Get results
    double res_cx,res_cy, res_cz,res_r;
    sphereFitting->getResult(res_cx, res_cy, res_cz, res_r);
     
    ASSERT_NEAR(    cx, res_cx, 1.0e-3);
    ASSERT_NEAR(    cy, res_cy, 1.0e-3);
    ASSERT_NEAR(    cz, res_cz, 1.0e-3);
    ASSERT_NEAR(radius,  res_r, 1.0e-3);
}

TEST(SphereFitting, isDegenerate)
{
    // Define estimation problem
    SphereFittingProblem * sphereFitting = new SphereFittingProblem();

    //sample 1
    std::vector<double> x1 = {0,1,2,3};
    std::vector<double> y1 = {0,1,2,3};
    std::vector<double> z1 = {0,1,2,3};
    sphereFitting->setData(x1,y1,z1);
    ASSERT_TRUE(sphereFitting->isDegenerate({0,1,2,3}));
    
    //sample 2
    std::vector<double> x2 = {0,1,2,3};
    std::vector<double> y2 = {0,1.001,2};
    std::vector<double> z2 = {0,1,2,0.001};
    sphereFitting->setData(x2,y2,z2);
    ASSERT_TRUE(sphereFitting->isDegenerate({0,1,2,3}));

    //sample 3
    std::vector<double> x3 = {0,1,5,7};
    std::vector<double> y3 = {0,1,2,4};
    std::vector<double> z3 = {0,1,1,2};
    sphereFitting->setData(x3,y3,z3);
    ASSERT_TRUE(sphereFitting->isDegenerate({0,1,2,3}));

    //sample 4
    std::vector<double> x4 = {0,  1,  0,   0};
    std::vector<double> y4 = {0,  0,  1,   0};
    std::vector<double> z4 = {1,  0,  0,  -1};
    sphereFitting->setData(x4,y4,z4);
    ASSERT_TRUE(!sphereFitting->isDegenerate({0,1,2,3}));

    //sample 5
    std::vector<double> x5 = {0.0,  0.0,  1.0,   0.0};
    std::vector<double> y5 = {0.0,  1.0,  0.0,   0.0};
    std::vector<double> z5 = {1.0,  0.0,  0.0,  -1.0};
    sphereFitting->setData(x5,y5,z5);
    ASSERT_TRUE(!sphereFitting->isDegenerate({0,1,2,3}));
}