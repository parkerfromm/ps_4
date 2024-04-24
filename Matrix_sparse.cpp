#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>
#include <array>
#include "myclasses.h"
#include "Matrix_dense.cpp"
using namespace std;

// template <typename T>
// class Matrix_sparse
// {

// public:
//     struct CSR
//     {
//         vector<T> values;
//         vector<int> column_indx;
//         vector<int> row_pointers; // indx of values that starts row of matrix

//         // if 5x5 matrix that means we have 25 enteries.
//         // non zero values will be stored in the values vector
//         /
//     };

//     Matrix_sparse();
//     Matrix_sparse(int rows, int cols);
//     Matrix_sparse(const Matrix_sparse<T> &other_matrix);
//     Matrix_sparse(int num_nonzero, int rows, int cols);

//     const T &operator()(int row, int col) const;

//     // use this one as fill function.
//     T &operator()(int row, int col);

//     vector<T> right_multiply(vector<T> x);
//     vector<T> left_multiply(vector<T> x);

//     array<int, 2> size() const;

//     Matrix_dense<T> convert_to_dense();

// private:
//     int _rows;
//     int _cols;
//     struct CSR _matrix;
// };

template <class T>
Matrix_sparse<T>::Matrix_sparse()
    : _rows(0), _cols(0) {}

template <typename T>
Matrix_sparse<T>::Matrix_sparse(int rows, int cols)
    : _rows(rows), _cols(cols)
{
    _matrix.row_pointers.assign(rows + 1, 0);

    // _matrix = std::vector<std::vector<T>>(_rows, std::vector<T>(_cols));
}

template <typename T>
Matrix_sparse<T>::Matrix_sparse(int num_nonzero, int rows, int cols)
    : _rows(rows), _cols(cols)
{
    _matrix.values.resize(num_nonzero);
    _matrix.column_indx.resize(num_nonzero);
    _matrix.row_pointers.assign(rows, 0);
};

// !!! update so you can create them for testing
// needs a method that can take the column row index and form matrix

// template <typename T>
// Matrix_sparse<T>::Matrix_sparse(vector<T> values,vector<int> column_indx, vector<int> row_pointers)
//     : _rows(rows), _cols(cols)
// {
//     _matrix.values= values;
//     _matrix.column_indx= column_idx;
//     _matrix.row_pointers.assign(rows, 0);
// };

template <typename T>
Matrix_sparse<T>::Matrix_sparse(const Matrix_sparse<T> &other_matrix)
    : _rows(other_matrix._rows), _cols(other_matrix._cols), _matrix(other_matrix._matrix) {}

// this is for reading and accessing

// need to edit this
// template <typename T>
// const T &Matrix_sparse<T>::operator()(int row, int col) const
// {

//     if (row >= _rows || col >= _cols || row < 0 || col < 0)
//     {
//         throw out_of_range("Index out of bounds");
//     }

//     int start_index = row_pointers[row];
//     int end_index = row_pointers[row + 1];

//     for (int j = start_index; j < end_indx; j++)
//     {

//         // add out of bound checks
//         if _matrix
//             .column_indx[j] == col
//             {
//                 return _matrix.values[j];
//             }
//         else
//         {
//             continue;
//         }
//     }

//     return 0;
// }

template <typename T>
const T &Matrix_sparse<T>::operator()(int row, int col) const
{
    if (row >= _rows || col >= _cols || row < 0 || col < 0)
    {
        throw std::out_of_range("Index out of bounds");
    }
    int start_index = _matrix.row_pointers[row];
    int end_index = _matrix.row_pointers[row + 1];

    for (int j = start_index; j < end_index; j++)
    {
        if (_matrix.column_indx[j] == col)
        {
            return _matrix.values[j];
        }
    }
    static T zero = T(); // Cache a static zero value to return
    return zero;
}

// to set A(2,3)=10.0
// template <typename T>
// T &Matrix_sparse<T>::operator()(int row, int col)
// {

//     // make sure number of rows and columns get declare matrix needs to be sized first;
//     // T c_ij;

//     if (row >= _rows || col >= _cols || row < 0 || col < 0)
//     {
//         throw out_of_range("Index out of bounds");
//     }

//     int start_index = row_pointers[row];
//     int end_index = row_pointers[row + 1];

//     for (int j = start_index; j < end_indx; j++)
//     {
//         if (_matrix.column_indx[i] == col)
//         {
//             return _matrix_values[j];
//         }
//         else if (_matrix.column_indx[i] > col)
//         {
//             _matrix.values.insert(_matrix.values.begin() + j, T{});
//             _matrix.column_indx.insert(_matrix.column_indx.begin() + j, col);

//             for (int i = row + 1; i < _rows + 1; i++)
//             {
//                 _matrix.row_pointers[i]++;
//             }

//             return _matrix_values[j];
//         }
//         // add out of bound checks
//         // check if need to fix the other version
//     }
// }

template <typename T>
T &Matrix_sparse<T>::operator()(int row, int col)
{
    if (row >= _rows || col >= _cols || row < 0 || col < 0)
    {
        throw std::out_of_range("Index out of bounds");
    }

    int start_index = _matrix.row_pointers[row];
    int end_index = _matrix.row_pointers[row + 1];

    for (int j = start_index; j < end_index; j++)
    {
        if (_matrix.column_indx[j] == col)
        {
            return _matrix.values[j];
        }
        else if (_matrix.column_indx[j] > col)
        {
            _matrix.values.insert(_matrix.values.begin() + j, T{});
            _matrix.column_indx.insert(_matrix.column_indx.begin() + j, col);
            // Update subsequent row pointers
            for (int i = row + 1; i <= _rows; i++)
            {
                _matrix.row_pointers[i]++;
            }
            return _matrix.values[j];
        }
    }

    // Handle case where element needs to be added at the end of the section
    _matrix.values.push_back(T{});
    _matrix.column_indx.push_back(col);
    _matrix.row_pointers[row + 1]++;
    for (int i = row + 2; i <= _rows; i++)
    {
        _matrix.row_pointers[i]++;
    }
    return _matrix.values.back();
}

// template <typename T>
// vector<T> Matrix_sparse<T>::left_multiply(vector<T> x)
// {

//     // mxn times nx1 gives mx1
//     vector<T> b(_rows);
//     for (int i = 0; i < _rows; i++)
//     {
//         int start_index = _matrix.row_pointers[i];
//         int end_index = _matrix.row_pointers[i + 1];

//         for (int j = start_index; j < end_indx; j++)
//         {
//             b[i] += _matrix.values[j] * x[_matrix.column_indx[i]];
//         }
//     }
// }

template <typename T>
vector<T> Matrix_sparse<T>::left_multiply(const vector<T> x) const
{
    vector<T> b(_rows, T{});
    for (int i = 0; i < _rows; i++)
    {
        int start_index = _matrix.row_pointers[i];
        int end_index = _matrix.row_pointers[i + 1];
        for (int j = start_index; j < end_index; j++)
        {
            b[i] += _matrix.values[j] * x[_matrix.column_indx[j]];
        }
    }
    return b;
}

// template <typename T>
// vector<T> Matrix_sparse<T>::right_multiply(vector<T> x)
// {

//     // 1x n nXm gives 1xm

//     vector<T> b(_cols);

//     for (int i = 0; i < _cols; i++)
//     {

//         for (int i = 0; i < _rows; i++)
//         {
//             int start_index = _matrix.row_pointers[i];
//             int end_index = _matrix.row_pointers[i + 1];

//             for (int k = start_index; k < end_indx; k++)
//             {
//                 b[_matrix.column_indx] += _matrix.values[k] * x[i];
//             }
//         }
//     }
// }

template <typename T>
vector<T> Matrix_sparse<T>::right_multiply(const vector<T> x) const
{
    vector<T> b(_cols, T{});
    for (int i = 0; i < _rows; i++)
    {
        int start_index = _matrix.row_pointers[i];
        int end_index = _matrix.row_pointers[i + 1];
        for (int j = start_index; j < end_index; j++)
        {
            b[_matrix.column_indx[j]] += _matrix.values[j] * x[i];
        }
    }
    return b;
}

template <typename T>
Matrix_dense<T> Matrix_sparse<T>::convert_to_dense()
{

    Matrix_dense<T> A(_rows, _cols);
    A.fill(0);

    for (int i = 0; i < _rows; i++)
    {
        int start_index = _matrix.row_pointers[i];
        int end_index = _matrix.row_pointers[i + 1];

        for (int j = start_index; j < end_index; j++)
        {
            A(i, _matrix.column_indx[j]) = _matrix.values[j];
        }
    }
    return A;
}

template <typename T>
vector<T> Matrix_sparse<T>::get_diag()
{

    vector<T> diags(_rows);
    // set this first equal to alll zeros then should be fine;
    // make it a class return

    for (int i = 0; i < _rows; i++)
    {
        int start_index = _matrix.row_pointers[i];
        int end_index = _matrix.row_pointers[i + 1];

        for (int j = start_index; j < end_index; j++)
        {
            if (_matrix.column_indx[j] == i)
            {

                diags[i] = _matrix.values[j];
            }
        }
    }

    return diags;
}
template <typename T>
int Matrix_sparse<T>::nrows() const
{
    return _rows;
}

template <typename T>
int Matrix_sparse<T>::ncols() const
{
    return _cols;
}
// template <typename T>
// array<int, 2> Matrix_sparse<T>::size() const
// {
//     return {_rows, _cols};
// }
template class Matrix_sparse<double>;
