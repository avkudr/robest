/**
 *  @brief Circle fitting to a set of 2D points
 *
 *  @author  Andrey Kudryavtsev (avkudr.github.io)
 *  @author  Rahima Djahel (github:rahma24000)
 *  @date    20/03/2018
 *  @version 1.0
 */

#include <list>
#include <vector>
#include <cmath>
#include <iostream>

#include "robust_estim.hpp"

struct Point2d{
    double x;
    double y;
};

typedef std::vector<Point2d> Point2Dvector;

class CircleFittingProblem : public robest::EstimationProblem{

public:
    CircleFittingProblem();
    ~CircleFittingProblem();

    void setData(std::vector<double> & x, std::vector<double> & y);

    double estimErrorForSample(int i);
    void   estimModelFromSamples(const std::vector<int> & samplesIdx);

    int getTotalNbSamples() const{
        return (int) points.size();
    }

    void getResult(double & res_cx, double & res_cy, double & res_r) const{
        res_cx = this->cx;
        res_cy = this->cy;
        res_r  = this->r;
    }

    inline bool isDegenerate(const std::vector<int> & samplesIdx)
    {
        Point2d & P = points[samplesIdx[0]];
        Point2d & V = points[samplesIdx[1]];
        Point2d & K = points[samplesIdx[2]];

        // verify that points P, V and K are not at the line -> verify that PV and PK are colinear:

        //1. calculate the directing coefficient of the line PV
        double f = (V.y-P.y)/(V.x-P.x);

        //2. calculate the directing coefficient of the line PK
        double h = (K.y-P.y)/(K.x-P.x);

        //3. PV and PK Are colineaire if and only if f = h
        return ( f - h < 1e-3 );
    }

private:

    Point2Dvector points; // input data
    double cx;
    double cy;
    double r; //radius
};
