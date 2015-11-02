#include <iostream>
#include <array>
#include <math.h>

#include <string.h>
#include <fstream>

using namespace std;

template <int dim>
class HyperBall
{
public:

    HyperBall(){};

    array<double,dim> evaluate(array<double,dim> & point);

private:
};
