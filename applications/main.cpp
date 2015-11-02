#include "matrix.h"
#include <typeinfo>

#include <boost/timer.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;
using namespace boost::numeric::ublas;

template<class T>
class BoostMatrix
{
public:
    BoostMatrix(unsigned int m,unsigned int n) : mat(m,n)
    {};

    BoostMatrix(const matrix<T> &in) : mat(in)
    {};

    void non_zero_init() {

        for (unsigned int i=0; i<mat.size1(); i++)
            for (unsigned int j=0; j<mat.size2(); j++)
                mat(i,j)=i*mat.size1()+j+1.;
    };

    void print()
    {
        cout << mat << endl;
    };

    matrix<T> get_boost_matrix() const
    {
        return mat;
    };

    const BoostMatrix<T> operator*(const BoostMatrix<T>& B) const
    {
        const auto b = B.get_boost_matrix();
        auto c = prod(mat, b);

        return BoostMatrix<T>(c);
    };


private:
    matrix<T> mat;
};

template<class M>
double run_test(const int m, const int k, const int l)
{
    M a(m,k);
    a.non_zero_init();

    M b(k,l);
    b.non_zero_init();

    M c(m,l);

    boost::timer t;
    c = a*b;
    return t.elapsed();

};

int main()
{

    cout << "boost matrix time = " << run_test<BoostMatrix<float>>(800,600,500) << endl;
    cout << "in house matrix time = " << run_test<Matrix<float>>(800,600,500) << endl;

    return 0;
}
