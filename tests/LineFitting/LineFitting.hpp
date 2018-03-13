#include <cmath>

#include<list>
#include<vector>
#include <cmath>

#include<list>
#include<vector>
#include<iostream>

#include"robust_estim.hpp"
using namespace std;


// ----------------------------------------------------------------------------
// ------- Class defining the Problem for RobustEstimator (RANSAC)
// ----------------------------------------------------------------------------
struct Point2d{
    double x;
    double y;
};

typedef std::vector<Point2d> Point2Dvector;

class LineFittingProblem : public EstimationProblem{

public:
    LineFittingProblem();
    ~LineFittingProblem();

    void setData(std::vector<double> & x, std::vector<double> & y);

    double estimErrorForSample(int i);
    double estimModelFromSamples(std::vector<int> samplesIdx);

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
