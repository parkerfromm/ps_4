
#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>
#include <array>

using namespace std;

#ifndef MYCLASSES_H
#define MYCLASSES_H

template <typename T>
class Matrix_dense
{
public:
    Matrix_dense();
    Matrix_dense(int rows, int cols);
    Matrix_dense(const Matrix_dense<T> &other_matrix);

    void resize(int rows, int cols);

    T &operator()(int row, int col);
    const T &operator()(int row, int col) const;

    Matrix_dense<T> multiply_by(const Matrix_dense<T> &other_matrix) const;

    Matrix_dense<T> operator+(const Matrix_dense<T> &other_matrix) const;
    Matrix_dense<T> operator*(const Matrix_dense<T> &other_matrix) const;

    Matrix_dense<T> &operator=(const Matrix_dense<T> &other_matrix);
    void fill(const T &value);
    Matrix_dense<T> &operator=(const T &value);
    bool operator==(const Matrix_dense<T> &other_matrix) const;

    template <typename t>
    friend ostream &operator<<(ostream &output, const Matrix_dense<t> &matrix);

    Matrix_dense<T> get_transpose();

    void transpose();

    std::array<int, 2> size() const;

    int nrows() const;
    int ncols() const;

    T *data();
    const T *data() const;

    void resize();
    void append_column(const std::vector<T> &new_column);

private:
    int _rows;
    int _cols;
    vector<vector<T>> _matrix;
    mutable std::vector<T> _contiguous_data;
};

template <typename T>
class Matrix_sparse
{

public:
    struct CSR
    {
        vector<T> values;
        vector<T> column_indx;
        vector<T> row_pointers; // indx of values that starts row of matrix

        // if 5x5 matrix that means we have 25 enteries.
        // non zero values will be stored in the values vector
    };

    Matrix_sparse();
    Matrix_sparse(int rows, int cols);
    Matrix_sparse(const Matrix_sparse<T> &other_matrix);
    Matrix_sparse(int num_nonzero, int rows, int cols);

    const T &operator()(int row, int col) const;

    // use this one as fill function.
    T &operator()(int row, int col);

    vector<T> right_multiply(const vector<T> x) const;
    vector<T> left_multiply(const vector<T> x) const;

    array<int, 2> size() const;

    Matrix_dense<T> convert_to_dense();
    vector<T> get_diag();
    int nrows() const;
    int ncols() const;

private:
    int _rows;
    int _cols;
    struct CSR _matrix;
};

#endif
// MYCLASSES_H