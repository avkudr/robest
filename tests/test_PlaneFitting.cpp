

#include "gtest/gtest.h"

#include <random>

#include "PlaneFitting/PlaneFitting.hpp"

void generateData(
    const std::vector<double> & modelParams,
    std::vector<double> & x, 
    std::vector<double> & y,
    std::vector<double> & z,
    double noiseVar = 0.0)
{

    double a = modelParams[0];
    double b = modelParams[1];
    double c = modelParams[2];
    double d = modelParams[3];

    double norm=sqrt(a*a+b*b+c*c+d*d);
    a = a / norm;
    b = b / norm;
    c = c / norm;
    d = d / norm;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0,noiseVar);

    for(int i = 0; i < 10; i++){
        for(int j = 0; j < 10; j++){
            x.push_back( i );
            y.push_back( j );
            z.push_back(( - a * x.back() - b * y.back() - d)/c);
        }
    }
}

TEST(PlaneFitting, idealCase)
{
    double a =  0.372997;
    double b = -0.136612;
    double c =  0.265316;
    double d = -0.878531;

    std::vector<double> x, y, z;
    generateData({a,b,c,d},x,y,z);

    PlaneFittingProblem * planeFitting = new PlaneFittingProblem();
    planeFitting->setData(x,y,z);

    //RANSAC * solver = new RANSAC();
    robest::MSAC * solver = new robest::MSAC();
    solver->solve(planeFitting);
    double res_a,res_b,res_c, res_d;
    planeFitting->getResult(res_a,res_b,res_c,res_d);

    //show result
    ASSERT_NEAR( fabs(a), fabs(res_a), 1.0e-6);
    ASSERT_NEAR( fabs(b), fabs(res_b), 1.0e-6);
    ASSERT_NEAR( fabs(c), fabs(res_c), 1.0e-6);
    ASSERT_NEAR( fabs(d), fabs(res_d), 1.0e-6);
}

TEST(PlaneFitting, idealCase2)
{
    double a = 0;
    double b = 1;
    double c = 5;
    double d = 0.58458;

    std::vector<double> x, y, z;
    generateData({a,b,c,d},x,y,z);

    PlaneFittingProblem * planeFitting = new PlaneFittingProblem();
    planeFitting->setData(x,y,z);

    //RANSAC * solver = new RANSAC();
    robest::MSAC * solver = new robest::MSAC();
    solver->solve(planeFitting);
    double res_a,res_b,res_c, res_d;
    planeFitting->getResult(res_a,res_b,res_c,res_d);

    double norm = sqrt(a*a+b*b+c*c+d*d);
    a = a / norm;
    b = b / norm;
    c = c / norm;
    d = d / norm;

    //show result
    ASSERT_NEAR( fabs(a), fabs(res_a), 1.0e-6);
    ASSERT_NEAR( fabs(b), fabs(res_b), 1.0e-6);
    ASSERT_NEAR( fabs(c), fabs(res_c), 1.0e-6);
    ASSERT_NEAR( fabs(d), fabs(res_d), 1.0e-6);
}   

TEST(PlaneFitting, idealCase3)
{
    double a = 0;
    double b = 0;
    double c = 1;
    double d = 0;

    std::vector<double> x, y, z;
    generateData({a,b,c,d},x,y,z);

    PlaneFittingProblem * planeFitting = new PlaneFittingProblem();
    planeFitting->setData(x,y,z);

    //RANSAC * solver = new RANSAC();
    robest::MSAC * solver = new robest::MSAC();
    solver->solve(planeFitting);
    double res_a,res_b,res_c, res_d;
    planeFitting->getResult(res_a,res_b,res_c,res_d);

    double norm = sqrt(a*a+b*b+c*c+d*d);
    a = a / norm;
    b = b / norm;
    c = c / norm;
    d = d / norm;

    //show result
    ASSERT_NEAR( a,      res_a, 1.0e-6);
    ASSERT_NEAR( b,      res_b, 1.0e-6);
    ASSERT_NEAR( c, fabs(res_c), 1.0e-6);
    ASSERT_NEAR( d,      res_d, 1.0e-6);
}
