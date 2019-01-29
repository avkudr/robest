/**
 *  @brief Line fitting to a set of 2D points
 *
 *  @author  Andrey Kudryavtsev (avkudr.github.io)
 *  @author  Rahima Djahel (github:rahma24000)
 *  @date    13/03/2018
 *  @version 1.0
 */

#include <cmath>
#include <list>
#include <vector>
#include <iostream>

#include "robust_estim.hpp"

struct Point2d{
    double x;
    double y;
};

typedef std::vector<Point2d> Point2Dvector;

class LineFittingProblem : public robest::EstimationProblem{

public:
    LineFittingProblem();
    ~LineFittingProblem();

    void setData(std::vector<double> & x, std::vector<double> & y);

    double estimErrorForSample(int i);
    void estimModelFromSamples(const std::vector<int> & samplesIdx);

    int getTotalNbSamples() const{
        return (int) points.size();
    }

    void getResult(double & resa, double & resb){
        resa = this->a;
        resb = this->b;
    }

private:
    Point2Dvector points; // Data
    double a;
    double b;

};
