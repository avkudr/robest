
#include "PlaneFitting.hpp"

PlaneFittingProblem::PlaneFittingProblem(){
    setNbParams(4);
    setNbMinSamples(3);
}

PlaneFittingProblem::~PlaneFittingProblem(){

}

void PlaneFittingProblem::setData(std::vector<double> & x, std::vector<double> & y,std::vector<double> & z){

    points.clear();
    for (auto i = 0; i < x.size(); i++){
        points.push_back(Point3d{x[i],y[i],z[i]});
    }
}

double PlaneFittingProblem::estimErrorForSample(int i)
{
    const Point3d & P = points[i];
    return std::fabs(A*P.x+B*P.y+C*P.z+D)/sqrt(A*A+B*B+C*C);
}

void PlaneFittingProblem::estimModelFromSamples(const std::vector<int> & samplesIdx){

    if( !isDegenerate(samplesIdx)){

        const Point3d & P = points[samplesIdx[0]];
        const Point3d & V = points[samplesIdx[1]];
        const Point3d & K = points[samplesIdx[2]];

        A =  (V.y-P.y)*(K.z-P.z) - (K.y-P.y)*(V.z-P.z);
        B =  (V.x-P.x)*(K.z-P.z) + (V.z-P.z)*(K.x-P.x);
        C =  (V.x-P.x)*(K.y-P.y) - (V.y-P.y)*(K.x-P.x);
        D =  (-1)*(A*P.x + B*P.y +C*P.z);
    }
}
bool PlaneFittingProblem::isDegenerate(std::vector<int> samplesIdx)
{
    const Point3d & P = points[samplesIdx[0]];
    const Point3d & V = points[samplesIdx[1]];
    const Point3d & K = points[samplesIdx[2]];

    const Point3d u{V.x - P.x, V.y - P.y, V.z - P.z};
    const Point3d v{K.x - P.x, K.y - P.y, K.z - P.z};

    const auto dotProduct = u.x * v.x + u.y * v.y + u.z * v.z;

    return ( dotProduct < 1e-5 );
}























