#include "LineFitting.hpp"

LineFittingProblem::LineFittingProblem()
{
    setNbParams(2);
    setNbMinSamples(2);
}

LineFittingProblem::~LineFittingProblem()
{

}

void LineFittingProblem::setData(std::vector<double> & x, std::vector<double> & y)
{
    points.clear();
    for (int i = 0; i < x.size(); i++){
        Point2d point;
        point.x=x[i];
        point.y=y[i];
        points.push_back(point);
    }
}

inline double LineFittingProblem::estimErrorForSample(int i)
{
    //distance point to line
    const Point2d & p = points[i];
    return std::fabs(a*p.x - p.y + b) / sqrt(a*a + 1); // distance line-point
}

inline void LineFittingProblem::estimModelFromSamples(const std::vector<int> & samplesIdx)
{
    //line from two points
    Point2d & P = points[samplesIdx[0]];
    Point2d & V = points[samplesIdx[1]];

    a = (V.y-P.y)/(V.x-P.x);
    b = P.y -a * P.x;
}




















