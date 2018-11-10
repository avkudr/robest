/**
 *  @brief Sphere fitting to a set of 3D points
 *
 *  @author  Andrey Kudryavtsev (avkudr.github.io)
 *  @author  Mark Anisimov (github:qM4RCp)
 *  @date    10/11/2018
 *  @version 1.0
 */

#include <cmath>
#include <list>
#include <vector>
#include <iostream>

#include "robust_estim.hpp"

struct Point3d{
    double x;
    double y;
    double z;
};

typedef std::vector<Point3d> Point3Dvector;

class SphereFittingProblem : public robest::EstimationProblem{

public:
    SphereFittingProblem();
    ~SphereFittingProblem();

    void setData(std::vector<double> & x, std::vector<double> & y, std::vector<double> & z);

    double estimErrorForSample(int i);
    void   estimModelFromSamples(std::vector<int> samplesIdx);

    double determinant(std::vector<int> samplesIdx);

    int getTotalNbSamples() const{
        return (int) points.size();
    }

    void getResult(double & res_cx, double & res_cy, double & res_cz, double & res_r) const{
        res_cx = this->cx;
        res_cy = this->cy;
        res_cz = this->cz;
        res_r  = this->r;
    }

private:
    Point3Dvector points;
    double cx;
    double cy;
    double cz;
    double r;
};