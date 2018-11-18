
#include "PlaneFitting.hpp"

PlaneFittingProblem::PlaneFittingProblem(){
    setNbParams(4);
    setNbMinSamples(3);
}

PlaneFittingProblem::~PlaneFittingProblem(){

}

void PlaneFittingProblem::setData(std::vector<double> & x, std::vector<double> & y,std::vector<double> & z){

    Point3d data;

    points.clear();

    for (int i = 0; i < x.size(); i++){
        data.x=x[i];
        data.y=y[i];
        data.z=z[i];

        points.push_back(data);
    }
}

double PlaneFittingProblem::estimErrorForSample(int i)
{
    Point3d & P = points[i];
    return std::abs(A*P.x+B*P.y+C*P.z+D)/sqrt(A*A+B*B+C*C);
}

void PlaneFittingProblem::estimModelFromSamples(std::vector<int> samplesIdx){

    if( !isDegenerate(samplesIdx)){

        Point3d & P = points[samplesIdx[0]];
        Point3d & V = points[samplesIdx[1]];
        Point3d & K = points[samplesIdx[2]];

        A =(V.y-P.y)*(K.z-P.z) - (K.y-P.y)*(V.z-P.z);
        B =-((V.x-P.x)*(K.z-P.z) - (V.z-P.z)*(K.x-P.x));
        C =(V.x-P.x)*(K.y-P.y) - (V.y-P.y)*(K.x-P.x);
        D =(-1)*(A*P.x + B*P.y +C*P.z);
    }
}
bool  PlaneFittingProblem::isDegenerate(std::vector<int> samplesIdx)
{
    Point3d & P = points[samplesIdx[0]];
    Point3d & V = points[samplesIdx[1]];
    Point3d & K = points[samplesIdx[2]];

    // verify that points P, V and K are not collinear ⇒ verify the angle between PV and PK == 0

    //F=PV.PK
    double F=(V.x-P.x * K.x-P.x)+(V.y-P.y * K.y-P.y)+(V.z-P.z * K.z-P.z);

    //H=||PV||
    double H=sqrt(std::abs(V.x-P.x * V.x-P.x)+(V.y-P.y * V.y-P.y)+(V.z-P.z * V.z-P.z));

    // L=||PK||
    double L=sqrt(std::abs(K.x-P.x * K.x-P.x)+(K.y-P.y * K.y-P.y)+(K.z-P.z * K.z-P.z));

    double theta=std::acos(F/H*L);

    return(theta<0.001);
}























