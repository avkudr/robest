#include "LineFitting.hpp"

using namespace std;

LineFittingProblem::LineFittingProblem(){
    setNbParams(2);
    setNbMinSamples(2);
}

LineFittingProblem::~LineFittingProblem(){

}

//
void LineFittingProblem::setData(std::vector<double> & x, std::vector<double> & y){

    Point2d data;

    for (int i = 0; i < x.size(); i++){
        data.x=x[i];
        data.y=y[i];

        points.push_back(data);
    }
}

double LineFittingProblem::estimErrorForSample(int i)
{
    Point2d & p = points[i];
    return std::fabs(a*p.x - p.y + b) / sqrt(a*a + 1); // distance line-point
}

double LineFittingProblem::estimModelFromSamples(std::vector<int> samplesIdx){

    Point2d & P = points[samplesIdx[0]];
    Point2d & V = points[samplesIdx[1]];

    a = (V.y-P.y)/(V.x-P.x);
    b = P.y -a * P.x;

    return 0;
}




















