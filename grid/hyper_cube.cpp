#include "hyper_cube.h"

template<int dim>
HyperCube<dim>::HyperCube(array<double,dim> & l, array<int,dim> & ni)
{
    length = l;
    n_intervals = ni;
    for (int i = 0; i < dim; i++)
    {
        double delta = length[i]/(double)n_intervals[i];

        for (int j =0; j < n_intervals[i]; j++)
            grid[i].push_back((double)j*delta);

        grid[i].push_back(length[i]);
    }
};

template<int dim>
void HyperCube<dim>::print()
{
    for (auto &row : grid)
        for (auto &kit : row)
            cout << kit << endl;
};


template<int dim>
void HyperCube<dim>::print_to_disk(string & grid_name)
{
    for (int i = 0; i < dim; i++)
    {
        string filename = grid_name + "_dim_" + to_string(i);
        auto row = grid[i];
        ofstream fout(filename, ios::out | ios::binary);
        fout.write((char*)&row[0], row.size() * sizeof(double));
        fout.close();
    }
};

template<int dim>
array<int,dim> HyperCube<dim>::flat_to_dim_id(int n_el)
{
    array<int,dim> dim_id;

    for (int i = dim; i > 1; i--)
    {
        int p = 1;
        for (int k=0; k<i-1; k++)
            p *= n_intervals[k];
        dim_id[i-1] = floor(n_el/p);
        n_el -= dim_id[i-1]*p;
    }
    dim_id[0] = n_el;

    return dim_id;
};

template class HyperCube<3>;
template class HyperCube<2>;
