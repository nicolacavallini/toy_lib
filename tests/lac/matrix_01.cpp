#include "matrix.h"
#include "../test.h"


using namespace std;

template<class T>
void run_test()
{
    Matrix<T> a(4,3);
    a.non_zero_init();
    output_file << "matrix a:" << endl;
    a.print_to_file(output_file);
    Matrix<T> b(3,5);
    b.non_zero_init();
    output_file << "matrix b:" << endl;
    b.print_to_file(output_file);
    
    auto c = a*b;

    output_file << "matrix c = a*b" << endl;
    c.print_to_file(output_file);
};

int main()
{
    run_test<double>();
    run_test<float>();
    return 0;
}
