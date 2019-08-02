#include "CircleFitting.hpp"

CircleFittingProblem::CircleFittingProblem(){
    setNbParams(3);
    setNbMinSamples(3);
}

CircleFittingProblem::~CircleFittingProblem(){

}

void CircleFittingProblem::setData(std::vector<double> & x, std::vector<double> & y)
{
    points.clear();
    for (int i = 0; i < x.size(); i++){
        Point2d data;
        data.x=x[i];
        data.y=y[i];
        points.push_back(data);
    }
}

inline double CircleFittingProblem::estimErrorForSample(int i)
{
    // distance circle-point = abs(<distance point-center> - radius)
    Point2d & p = points[i];
    return std::abs(sqrt((p.x-cx)*(p.x-cx)+(p.y-cy)*(p.y-cy)) - r);
}

inline void CircleFittingProblem::estimModelFromSamples(const std::vector<int> & samplesIdx){
    if( !isDegenerate(samplesIdx)){
        Point2d & P = points[samplesIdx[0]];
        Point2d & V = points[samplesIdx[1]];
        Point2d & K = points[samplesIdx[2]];
        //calculation of the coefficients of the mediating straight lines
        double a = -(V.x - P.x)/(V.y - P.y);
        double b = (V.x * V.x - P.x * P.x + V.y * V.y - P.y * P.y)/(2* (V.y - P.y));
        double c = -(K.x - V.x)/(K.y - V.y);
        double d = (K.x * K.x - V.x * V.x + K.y * K.y - V.y * V.y)/(2* (K.y - V.y));

        //calculate the coordinates of the center of the circle O(A,B)
        cx = (b-d)/(c-a);
        cy = a*cx + b;

        //calculate the radius of a circle
        r = sqrt((P.x - cx)*(P.x - cx)+(P.y - cy)*(P.y - cy));
    }
}
























