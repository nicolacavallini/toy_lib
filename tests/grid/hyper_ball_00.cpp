#include <math.h>
#include <array>


#include "hyper_ball.h"
#include "../test.h"


using namespace std;

template<int dim>
void run_test(array<double,dim> &p)
{

    HyperBall<dim> h_ball;
    
    auto m = h_ball.evaluate(p);
       
    output_file << "testing point: ";
    for (auto &c : m)
            output_file << c << ", ";
    output_file << endl;
    
    output_file << "resulting mapped point: ";
    for (auto &c : m)
            output_file << c << ", ";
    output_file << endl;
}

int main()
{
    double t = M_PI/4;
    double r = 2/sqrt(2);
    array<double,2> p = {r,t};
    
    run_test<2>(p);
    
    double f = M_PI/4;
    array<double,3> q = {r,t,f};
    
    run_test<3>(q);
    

 
    return 0;
}
