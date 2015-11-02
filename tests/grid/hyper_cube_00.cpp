#include <iostream>
#include <vector>
#include <array>
#include <math.h>

#include <string.h>
#include <fstream>

#include "hyper_cube.h"
#include "../test.h"


using namespace std;

template<int dim>
void run_test(array<double,dim> &length, array<int,dim> & n_intervals, int element)
{
    HyperCube<dim> h_cube(length,n_intervals);
    
    string log =  "dim = " + to_string(dim) + ", testing element = " + to_string(element);
    
    output_file << log <<endl;

    auto dim_id = h_cube.flat_to_dim_id(element);

    for (int i = 0; i < dim; i++)
        output_file << "id[" << i << "] = "<< dim_id[i] <<endl;
    
    output_file << endl;
}

int main()
{
    array<double,3> length = {1,2,3};
    array<int,3> n_intervals = {2,4,3};
    
    run_test<3>(length,n_intervals,16);
    run_test<3>(length,n_intervals,21);
    
    array<double,2> l = {2,3};
    array<int,2> ni = {4,3};
    
    run_test<2>(l,ni,6);
    run_test<2>(l,ni,9);
        
    return 0;
}
