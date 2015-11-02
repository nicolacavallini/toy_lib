#ifndef __matrix_H
#define __matrix_H
#ifdef USEBLAS
#include <cblas.h>
#endif
#include <iostream>

template <class T>
class Matrix
{
	private:
		T* data;
		T dump;

	protected:
        int m_n, m_m;
		Matrix();
      
	public:
		virtual const T trace(void) const { return 0.0f; };
		
		Matrix(int, int);
		
		virtual ~Matrix();
		
		Matrix(const Matrix<T> & mat);
		
		int get_n(void) const;
		
		int get_m(void) const;
		
		void print(void) const;
		void print_to_file(std::ofstream &output_file) const;
        
        Matrix<T> &operator=(const Matrix<T>& );
		
		const Matrix<T> operator*(const Matrix<T>&) const;
		
		T* operator[](int ) const;
		
		const T &operator()(int, int ) const;
		
		T &operator()(int, int );
		
		void non_zero_init();
};


#endif
