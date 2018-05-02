/**
 *  @brief Line fitting to a set of 2D points
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
    void   estimModelFromSamples(std::vector<int> samplesIdx);

    int getTotalNbSamples() const{
        return (int) points.size();
    }

    void getResult(double & res_cx, double & res_cy, double & res_r) const{
        res_cx = this->cx;
        res_cy = this->cy;
        res_r  = this->r;
    }

    bool isDegenerate(std::vector<int> samplesIdx);

private:

    Point2Dvector points; // input data
    double cx;
    double cy;
    double r; //radius
};
