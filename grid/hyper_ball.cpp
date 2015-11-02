#include "hyper_ball.h"

template<>
array<double,2> HyperBall<2>::evaluate(array<double,2> & point)
{
    double r = point[0];
    double theta = point[1];

    array<double,2> map;

    map[0] = r*cos(theta);
    map[1] = r*sin(theta);

    return map;
}

template<>
array<double,3> HyperBall<3>::evaluate(array<double,3> & point)
{
    double r = point[0];
    double theta = point[1];
    double phi = point[2];

    array<double,3> map;

    map[0] = r*sin(theta)*cos(phi);
    map[1] = r*sin(theta)*sin(phi);
    map[2] = r*cos(theta);

    return map;
}

template class HyperBall<2>;
template class HyperBall<3>;
