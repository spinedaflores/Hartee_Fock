#ifndef MY_MATRIX_H
#define MY_MATRIX_H
#include <vector>
#include <iostream>
#include <iomanip>


template <typename T>
struct matrix
{
	int col, row;
	T* data;

	matrix(int c, int r) : col(c), row(r) 
	{
		data = new T[col*row]; 
		std::fill(data, data+(row*col), 0);
	}

	~matrix()
	{
		delete [] data;
	}


	T& operator () (const int i, const int j)
	{
		return data[i+(row*j)];
	}

	T& operator [] (const int i)
	{
		return data[i];
	}

	T* begin()
	{
		return &(data[0]);
	}

	T* end()
	{
		return &(data[col*row]);
	}

	int size()
	{
		return (row*col);
	}

	int rows()
	{
		return row;
	}

	int cols()
	{
		return col;
	}

	void print()
	{
		for (int i = 0; i < row; i++)
		{
				for (int j = 0; j < col; j++)
				{
						if (j != 0) 
							std::cout << " ";
						std::cout << std::setw(5) << std::scientific << std::setw(10) << std::right << std::setfill(' ') << std::setprecision(5) << data[i+(row*j)];
				}
				std::cout << "\n";
		}
	}

};



#endif
