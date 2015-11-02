#include <iostream>
#include <vector>
#include <array>
#include <math.h>

#include <string.h>
#include <fstream>

using namespace std;

template <int dim>
class HyperCube
{
public:

    HyperCube(array<double,dim> & l, array<int,dim> & ni);

    void print();
 
    void print_to_disk(string & grid_name);

    array<int,dim> flat_to_dim_id(int n_el);


private:
    array<vector<double>,dim> grid;
    array<double,dim> length;
    array<int,dim> n_intervals;
};
