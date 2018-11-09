#include "gtest/gtest.h"

#include <random>

#include "SphereFitting/SphereFitting.hpp"

//gaussian noise generation
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
    std::vector<double> x = {0,  2,  0,   0};
    std::vector<double> y = {0,  0,  2,   0};
    std::vector<double> z = {2,  0,  0,  -2};

    double cx = 0; // sphere center C:(cx,cy, cz)
    double cy = 0;
    double cz = 0;
    double radius = 2; //sphere radius
    double noiseVar = 0.0;

    SphereFittingProblem * sphereFitting = new SphereFittingProblem();
    sphereFitting->setData(x,y,z);

    robest::MSAC * MSACsolver = new robest::MSAC();
    MSACsolver->solve(sphereFitting);

    double res_cx,res_cy, res_cz,res_r;
    sphereFitting->getResult(res_cx, res_cy, res_cz, res_r);

    ASSERT_NEAR(    cx, res_cx, 1.0e-11);
    ASSERT_NEAR(    cy, res_cy, 1.0e-11);
    ASSERT_NEAR(    cz, res_cz, 1.0e-11);
    ASSERT_NEAR(radius,  res_r, 1.0e-11);
}   

TEST(SphereFitting, idealCase2)
{
    std::vector<double> x = {8,     12,     4,      8,      8,      8};
    std::vector<double> y = {10,    6,      6,      2,      6,      6};
    std::vector<double> z = {-11,   -11,    -11,    -11,    -15,    -7};

    double cx = 8; // sphere center C:(cx,cy,cz)
    double cy = 6;
    double cz = -11;
    double radius = 4; //sphere radius
    double noiseVar = 0.0;

    //generateSphereData(cx,cy,cz,radius,noiseVar,x,y,z);

    SphereFittingProblem * sphereFitting = new SphereFittingProblem();
    sphereFitting->setData(x,y,z);

    robest::MSAC * MSACsolver = new robest::MSAC();
    MSACsolver->solve(sphereFitting);

    double res_cx,res_cy, res_cz,res_r;
    sphereFitting->getResult(res_cx, res_cy, res_cz, res_r);
    /*
    std::cout << std::endl;
    std::cout << "Estimated centre of shpere: cx = " << res_cx << ", cy = " << res_cy << ", cz = " << res_cz << std::endl;
    std::cout << "Estimated radius of shpere: r = " << res_r  << std::endl;
    std::cout << std::endl;
    */
    ASSERT_NEAR(    cx, res_cx, 1.0e-11);
    ASSERT_NEAR(    cy, res_cy, 1.0e-11);
    ASSERT_NEAR(    cz, res_cz, 1.0e-11);
    ASSERT_NEAR(radius,  res_r, 1.0e-11);
}   

TEST(SphereFitting, smallNoise)
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    double cx = 5.23542; //sphere center C:(cx,cy,cz)
    double cy = 3.35662;
    double cz = -1.96731;
    double radius = 4.54854; //sphere radius
    double noiseVar = 0.001;

    generateSphereData(cx,cy,cz,radius,noiseVar,x,y,z);

    SphereFittingProblem * sphereFitting = new SphereFittingProblem();
    sphereFitting->setData(x,y,z);

    robest::MSAC * MSACsolver = new robest::MSAC();
    MSACsolver->solve(sphereFitting);

    double res_cx,res_cy, res_cz,res_r;
    sphereFitting->getResult(res_cx, res_cy, res_cz, res_r);
    /*
    std::cout << std::endl;
    std::cout << "Estimated centre of shpere: cx = " << res_cx << ", cy = " << res_cy << ", cz = " << res_cz << std::endl;
    std::cout << "Estimated radius of shpere: r = " << res_r  << std::endl;
    std::cout << std::endl;
    */
    ASSERT_NEAR(    cx, res_cx, 1.0e-3);
    ASSERT_NEAR(    cy, res_cy, 1.0e-3);
    ASSERT_NEAR(    cz, res_cz, 1.0e-3);
    ASSERT_NEAR(radius,  res_r, 1.0e-3);
}

TEST(SphereFitting, outliers)
{
    std::vector<double> x = {3, 8,  -2, 3,  3,  3,  16, -24, 0};
    std::vector<double> y = {3, -2, -2, -7, -2, -2, 12, 16, -15, 0};
    std::vector<double> z = {7, 7,  7,  7,  2,  12, -2, 3,  0};

    double cx = 3; // sphere center C:(cx,cy,cz)
    double cy = -2;
    double cz = 7;
    double radius = 5; //sphere radius;

    SphereFittingProblem * sphereFitting = new SphereFittingProblem();
    sphereFitting->setData(x,y,z);

    robest::MSAC * MSACsolver = new robest::MSAC();
    MSACsolver->solve(sphereFitting);

    double res_cx,res_cy, res_cz,res_r;
    sphereFitting->getResult(res_cx, res_cy, res_cz, res_r);
    /*
    std::cout << std::endl;
    std::cout << "Estimated centre of shpere: cx = " << res_cx << ", cy = " << res_cy << ", cz = " << res_cz << std::endl;
    std::cout << "Estimated radius of shpere: r = " << res_r  << std::endl;
    std::cout << std::endl;
    */
    ASSERT_NEAR(    cx, res_cx, 1.0e-3);
    ASSERT_NEAR(    cy, res_cy, 1.0e-3);
    ASSERT_NEAR(    cz, res_cz, 1.0e-3);
    ASSERT_NEAR(radius,  res_r, 1.0e-3);
}