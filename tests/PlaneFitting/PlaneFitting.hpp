/**
 *  @brief Plane fitting to a set of 3D points
 *
 *  @author  Rahima Djahel (github:rahma24000)
 *  @author  Andrey Kudryavtsev (avkudr.github.io)
 *  @date    18/11/2018
 *  @version 1.0
 */

#include <iostream>
#include <vector>

#include "robust_estim.hpp"

struct Point3d{
    double x;
    double y;
    double z;
};

typedef std::vector<Point3d> Point3Dvector;

class PlaneFittingProblem : public robest::EstimationProblem{

public:
    PlaneFittingProblem();
    ~PlaneFittingProblem();

    void setData(std::vector<double> & x, std::vector<double> & y,std::vector<double> & z );

    double estimErrorForSample(int i);
    void estimModelFromSamples(const std::vector<int> & samplesIdx);

    int getTotalNbSamples() const{
        return (int) points.size();
    }
    void getResult(double & resa, double & resb, double & resc, double &resd){

        const double & a = this->A;
        const double & b = this->B;
        const double & c = this->C;
        const double & d = this->D;
        double norm = sqrt(a*a+b*b+c*c+d*d);
        if (norm == 0) norm = 1;

        resa = a / norm;
        resb = b / norm;
        resc = c / norm;
        resd = d / norm;
    }
private:
    bool isDegenerate(std::vector<int> samplesIdx);
    Point3Dvector points; // Data
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double D = 0.0;
};
