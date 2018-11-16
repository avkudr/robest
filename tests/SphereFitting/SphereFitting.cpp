#include "SphereFitting.hpp"

SphereFittingProblem::SphereFittingProblem(){
    setNbParams(4);
    setNbMinSamples(4);
}

SphereFittingProblem::~SphereFittingProblem(){

}

void SphereFittingProblem::setData(std::vector<double> & x, std::vector<double> & y, std::vector<double> & z)
{
    points.clear();
    for (int i = 0; i < x.size(); i++){
        Point3d point;
        point.x=x[i];
        point.y=y[i];
        point.z=z[i];
        points.push_back(point);
    }
}

inline double SphereFittingProblem::estimErrorForSample(int i)
{
    const Point3d & p = this->points[i];
    return std::fabs(std::sqrt(std::pow((p.x - this->cx),2)+ std::pow((p.y - this->cy),2) + std::pow((this->cz - p.z),2)) - this->r);
}

inline void SphereFittingProblem::estimModelFromSamples(std::vector<int> samplesIdx)
{
    const Point3Dvector & P = this->points;

    const auto & x1 = P[samplesIdx[0]].x, x2 = P[samplesIdx[1]].x, x3 = P[samplesIdx[2]].x, x4 = P[samplesIdx[3]].x;
    const auto & y1 = P[samplesIdx[0]].y, y2 = P[samplesIdx[1]].y, y3 = P[samplesIdx[2]].y, y4 = P[samplesIdx[3]].y;
    const auto & z1 = P[samplesIdx[0]].z, z2 = P[samplesIdx[1]].z, z3 = P[samplesIdx[2]].z, z4 = P[samplesIdx[3]].z;

    double D = this->determinantFromDataPoints(samplesIdx);

    if (!(std::fabs(D) < 1e-3)) // equivalent of isDegenerate
    {    
        double f1 = std::pow(x1,2)+std::pow(y1,2)+std::pow(z1,2);
        double f2 = std::pow(x2,2)+std::pow(y2,2)+std::pow(z2,2);
        double f3 = std::pow(x3,2)+std::pow(y3,2)+std::pow(z3,2);
        double f4 = std::pow(x4,2)+std::pow(y4,2)+std::pow(z4,2);

        //calculation of the centre of the sphere
        this->cx = (f3*(4*(y1*(z2 - z4) + y2*(z4 - z1) + y4*(z1 - z2)))\
                    - f4*(4*(y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)))\
                    - f2*(4*(y1*(z3 - z4) + y3*(z4 - z1) + y4*(z1 - z3)))\
                    + f1*(4*(y2*(z3 - z4) + y3*(z4 - z2) + y4*(z2 - z3)))) / D;
       
        this->cy = (f4*(4*(x1*(z2 - z3) + x2*(z3 - z1) + x3*(z1 - z2)))\
                    - f3*(4*(x1*(z2 - z4) + x2*(z4 - z1) + x4*(z1 - z2)))\
                    + f2*(4*(x1*(z3 - z4) + x3*(z4 - z1) + x4*(z1 - z3)))\
                    - f1*(4*(x2*(z3 - z4) + x3*(z4 - z2) + x4*(z2 - z3)))) / D;
        
        this->cz = (f3*(4*(x1*(y2 - y4) + x2*(y4 - y1) + x4*(y1 - y2)))\
                    - f4*(4*(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)))\
                    - f2*(4*(x1*(y3 - y4) + x3*(y4 - y1) + x4*(y1 - y3)))\
                    + f1*(4*(x2*(y3 - y4) + x3*(y4 - y2) + x4*(y2 - y3)))) / D;
        
        //calculation of the radius of the sphere
        double co = (f4*(8*(x1*(y2*z3 - y3*z2) + x2*(y3*z1 - y1*z3) + x3*(y1*z2 - y2*z1)))\
                    - f3*(8*(x1*(y2*z4 - y4*z2) + x2*(y4*z1 - y1*z4) + x4*(y1*z2 - y2*z1)))\
                    + f2*(8*(x1*(y3*z4 - y4*z3) + x3*(y4*z1 - y1*z4) + x4*(y1*z3 - y3*z1)))\
                    - f1*(8*(x2*(y3*z4 - y4*z3) + x3*(y4*z2 - y2*z4) + x4*(y2*z3 - y3*z2)))) / D;
        
        this->r = std::sqrt(std::pow(this->cx,2) + std::pow(this->cy,2) + std::pow(this->cz,2) + co);
    } 
}

inline double SphereFittingProblem::determinantFromDataPoints(const std::vector<int> & samplesIdx)
{
    const Point3Dvector & P = this->points;

    const auto & x1 = P[samplesIdx[0]].x, x2 = P[samplesIdx[1]].x, x3 = P[samplesIdx[2]].x, x4 = P[samplesIdx[3]].x;
    const auto & y1 = P[samplesIdx[0]].y, y2 = P[samplesIdx[1]].y, y3 = P[samplesIdx[2]].y, y4 = P[samplesIdx[3]].y;
    const auto & z1 = P[samplesIdx[0]].z, z2 = P[samplesIdx[1]].z, z3 = P[samplesIdx[2]].z, z4 = P[samplesIdx[3]].z;
    
    return  8*(x1*(y2*(z3 - z4) + y3*(z4 - z2) + y4*(z2 - z3))\
             + x2*(y1*(z4 - z3) + y3*(z1 - z4) + y4*(z3 - z1))\
             + x3*(y1*(z2 - z4) + y2*(z4 - z1) + y4*(z1 - z2))\
             + x4*(y1*(z3 - z2) + y2*(z1 - z3) + y3*(z2 - z1)));
}

bool SphereFittingProblem::isDegenerate(const std::vector<int> & samplesIdx)
{
    return (std::fabs(this->determinantFromDataPoints(samplesIdx)) < 1e-3);
}