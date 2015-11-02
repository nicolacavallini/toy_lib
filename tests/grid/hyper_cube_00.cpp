#include <iostream>
#include <vector>
#include <array>
#include <math.h>

#include <string.h>
#include <fstream>

#include "hyper_cube.h"
#include "../test.h"


using namespace std;

int main()
{
    const int dim = 3;
    array<double,dim> length = {1,2,3};
    array<int,dim> n_intervals = {2,4,3};
    

    HyperCube<dim> h_cube(length,n_intervals);
    
    output_file << "testing element 16:" <<endl;

    auto dim_id = h_cube.flat_to_dim_id(16);

    for (int i = 0; i < dim; i++)
        output_file << "id[" << i << "] = "<< dim_id[i] <<endl;

    output_file << "testing element 21:" << endl;

    dim_id = h_cube.flat_to_dim_id(21);

    for (int i = 0; i < dim; i++)
        output_file << "id[" << i << "] = "<< dim_id[i] <<endl;

    return 0;
}
