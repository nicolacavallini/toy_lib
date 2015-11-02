#include "matrix.h"
#include <fstream>

template<class T>
Matrix< T > ::Matrix()
{
#ifdef DEBUG
        std::cout << " default Matrix constructor " << std::endl;
#endif
}

template<class T>
Matrix< T >::Matrix(int n, int m) 
:
m_n(n),
m_m(m)
{
#ifdef DEBUG
        std::cout << " custom Matrix constructor " << std::endl;
#endif
	data=new T[m_m*m_n];
}

template<class T>
Matrix< T >::~Matrix()
{
#ifdef DEBUG
        std::cout << " Matrix destructor " << std::endl;
#endif
	delete [] data;
}

template<class T>
Matrix<T>::Matrix(const Matrix<T> & mat)
{
#ifdef DEBUG
	std::cout << " Matrix copy constructor" << std::endl; 
#endif
	m_m=mat.get_m();
	m_n=mat.get_n();
        data=new T[m_m*m_n];
	for (int i=0; i<m_n;i++)
		for (int j=0; j<m_m; j++)
			data[i*m_m+j]=mat.data[i*m_m+j];
}

template<class T>
int Matrix<T>::get_n() const
{
  return m_n;
}


template<class T>
void Matrix<T>::non_zero_init()
{
	for (int i=0; i<m_n;i++)
		for (int j=0; j<m_m; j++)
			data[i*m_m+j]=i*m_m+j+1.;
}

template<class T>
int Matrix<T>::get_m() const
{
  return m_m;
}

template<class T>
void Matrix<T>::print_to_file(std::ofstream &output_file) const
{
	for (int i=0; i<m_n;i++)
		{
		for (int j=0; j<m_m; j++)
                	output_file << " " << data[i*m_m+j];
      		output_file<<std::endl;
		}
}

template<class T>
void Matrix<T>::print() const
{
	for (int i=0; i<m_n;i++)
		{
		for (int j=0; j<m_m; j++)
                	std::cout << " " << data[i*m_m+j];
      		std::cout<<std::endl;
		}
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix& mat)
{
	if (this != &mat)
		{
#ifdef DEBUG
		std::cout << " Matrix assignment operator " << std::endl;
#endif
		delete [] data;
		m_m=mat.m_m;
		m_n=mat.m_n;
		data=new T[m_m*m_n];
        	for (int i=0; i<m_n;i++)
                	for (int j=0; j<m_m; j++)
                        	data[i*m_m+j]=mat.data[i*m_m+j];
		}
	return *this;	
}

template<class T>
const Matrix<T> Matrix<T>::operator*(const Matrix<T>& mat) const
{
	Matrix c(get_n(), mat.get_m());
	if (get_m() != mat.get_n())
		{
		std::cerr <<"incompatible shapes" << std::endl;
		for(int i=0; i < c.get_n()*c.get_m();i++)
			c.data[i] = 0.0f;
		return c;
		}

	T temp;
	for(int i=0; i < get_n();i++)
		for (int j=0; j < mat.get_m(); j++)
			{
			temp = 0.0f;
			for (int k=0; k< get_m(); k++)
				temp += data[i*get_m()+k] * mat.data[k*mat.get_m()+j];
			c.data[i*c.get_m()+j] = temp;
			}
	return c;		
}

template<class T>
T* Matrix<T>::operator[](int i) const
{
	return &data[i*m_m];
}

template<class T>
const T& Matrix<T>::operator()(int i, int j) const 
{
#if DEBUG == 2 
	std::cout << " const Matrix operator()(int i, int j)  " << std::endl;
#endif
	if (i>=0 and j>=0 and i<m_n and j<m_m)
		{
		return data[i*m_m+j];
		}
	else
		{
		std::cerr << " index out of bound! " << std::endl;
		return dump; 
		}
}

template<class T>
T& Matrix<T>::operator()(int i, int j)
{
#if DEBUG == 2
	std::cout << " non-const Matrix operator()(int i, int j) " << std::endl;
#endif
        if (i>=0 and j>=0 and i<m_n and j<m_m)
                {
                return data[i*m_m+j];
                }
        else
                {
                std::cerr << " index out of bound! " << std::endl;
                return dump;
                }
}


template<>
const Matrix<double> Matrix<double>::operator*(const Matrix<double>& mat) const
{
        Matrix c(get_n(), mat.get_m());
        if (get_m() != mat.get_n())
                {
                std::cerr <<"incompatible shapes" << std::endl;
                for(int i=0; i < c.get_n()*c.get_m();i++)
                        c.data[i] = 0.0;
                return c;
                }

#ifdef USEBLAS
	std::cout << " specialized operator* - using cblas_dgemm " << std::endl;
	double alpha=1.0;
	double beta=0.0;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_n, mat.m_m, m_m, alpha, data, m_m, mat.data, mat.m_m, beta, c.data, c.m_m);
#else
	std::cout << " specialized operator* - using local implementation for double " << std::endl;
	double temp;
        for(int i=0; i < get_n();i++)
                for (int j=0; j < mat.get_m(); j++)
                        {
                        temp = 0.0;
                        for (int k=0; k< get_m(); k++)
                                temp += data[i*get_m()+k] * mat.data[k*mat.get_m()+j];
                        c.data[i*c.get_m()+j] = temp;
                        }
#endif
        return c;
}

template<>
const Matrix<float> Matrix<float>::operator*(const Matrix<float>& mat) const
{
        Matrix c(get_n(), mat.get_m());
        if (get_m() != mat.get_n())
                {
                std::cerr <<"incompatible shapes" << std::endl;
                for(int i=0; i < c.get_n()*c.get_m();i++)
                        c.data[i] = 0.0f;
                return c;
                }

#ifdef USEBLAS
        std::cout << " specialized operator* - using cblas_sgemm " << std::endl;
        float alpha=1.0f;
        float beta=0.0f;
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_n, mat.m_m, m_m, alpha, data, m_m, mat.data, mat.m_m, beta, c.data, c.m_m);
#else
        std::cout << " specialized operator* - using local implementation for float " << std::endl;
        float temp;
        for(int i=0; i < get_n();i++)
                for (int j=0; j < mat.get_m(); j++)
                        {
                        temp = 0.0f;
                        for (int k=0; k< get_m(); k++)
                                temp += data[i*get_m()+k] * mat.data[k*mat.get_m()+j];
                        c.data[i*c.get_m()+j] = temp;
                        }
#endif
        return c;
}

template class Matrix<double>;
template class Matrix<float>;
