#include <vector>

using namespace std;

#include <vector>
#include <iostream>
#include <complex>
#include <algorithm>
#include <array>
#include <numeric>
#include <random>
#include "myclasses.h"
using namespace std;

// template <typename T>
//  class Matrix_dense
//  {
//  public:
//      Matrix_dense();
//      Matrix_dense(int rows, int cols);
//      Matrix_dense(const Matrix_dense<T> &other_matrix);

//     void resize(int rows, int cols);

//     T &operator()(int row, int col);
//     const T &operator()(int row, int col) const;

//     Matrix_dense<T> multiply_by(const Matrix_dense<T> &other_matrix) const;

//     Matrix_dense<T> operator+(const Matrix_dense<T> &other_matrix) const;
//     Matrix_dense<T> operator*(const Matrix_dense<T> &other_matrix) const;

//     Matrix_dense<T> &operator=(const Matrix_dense<T> &other_matrix);
//     void fill(const T &value);
//     Matrix_dense<T> &operator=(const T &value);
//     bool operator==(const Matrix_dense<T> &other_matrix) const;

//     template <typename t>
//     friend ostream &operator<<(ostream &output, const Matrix_dense<t> &matrix);

//     Matrix_dense<T> get_transpose();

//     void transpose();

//     std::array<int, 2> size() const;

//     const int nrows() const;
//     const int ncols() const;

//     T *data();

// private:
//     int _rows;
//     int _cols;
//     vector<vector<T>> _matrix;
// };

template <typename T>
Matrix_dense<T>::Matrix_dense()
    : _rows(0), _cols(0) {}

template <typename T>
Matrix_dense<T>::Matrix_dense(int rows, int cols)
    : _rows(rows), _cols(cols)
{
    _matrix = std::vector<std::vector<T>>(_rows, std::vector<T>(_cols));
}

template <typename T>
Matrix_dense<T>::Matrix_dense(const Matrix_dense<T> &other_matrix)
    : _rows(other_matrix._rows), _cols(other_matrix._cols), _matrix(other_matrix._matrix) {}

// set individual componet to value ex. m(2,3) = 10.0
template <typename T>
T &Matrix_dense<T>::operator()(int row, int col)
{
    if (row >= _rows || col >= _cols || row < 0 || col < 0)
    {
        throw out_of_range("Index out of bounds");
    }
    return _matrix[row][col];
}

// access (read) the component
template <typename T>
const T &Matrix_dense<T>::operator()(int row, int col) const
{
    if (row >= _rows || col >= _cols || row < 0 || col < 0)
    {
        throw out_of_range("Index out of bounds");
    }
    return _matrix[row][col];
}

template <typename T>
void Matrix_dense<T>::resize(int rows, int cols)
{
    if (rows < 0 || cols < 0)
    {
        throw std::invalid_argument("Matrix_dense size cannot be negative.");
    }
    _rows = rows;
    _cols = cols;
    _matrix = std::vector<std::vector<T>>(_rows, std::vector<T>(_cols));
}

template <typename T>
Matrix_dense<T> Matrix_dense<T>::operator+(const Matrix_dense<T> &other_matrix) const
{
    if (_rows != other_matrix._rows || _cols != other_matrix._cols)
    {
        throw invalid_argument("Matrices dimensions must match for addition");
    }

    Matrix_dense<T> new_matrix(_rows, _cols);
    for (int i = 0; i < _rows; i++)
    {
        for (int j = 0; j < _cols; j++)
        {
            new_matrix._matrix[i][j] = _matrix[i][j] + other_matrix._matrix[i][j];
        }
    }
    return new_matrix;
}

template <typename T>
Matrix_dense<T> Matrix_dense<T>::operator*(const Matrix_dense<T> &other_matrix) const
{
    if (_cols != other_matrix._rows)
    {
        throw invalid_argument("Matrices columns and rows dimensions must match for multiplication");
    }

    Matrix_dense<T> new_matrix(_rows, other_matrix._cols);
    for (int i = 0; i < _rows; i++)
    {
        for (int j = 0; j < other_matrix._cols; j++)
        {
            for (int k = 0; k < _cols; k++)
            {
                new_matrix._matrix[i][j] += _matrix[i][k] * other_matrix._matrix[k][j];
            }
        }
    }
    return new_matrix;
}

template <typename T>
Matrix_dense<T> Matrix_dense<T>::multiply_by(const Matrix_dense<T> &other_matrix) const
{
    if (_rows != other_matrix._rows || _cols != other_matrix._cols)
    {
        throw std::invalid_argument("Matrices dimensions must match for element wise multiplication");
    }

    Matrix_dense<T> new_matrix(_rows, _cols);
    for (int i = 0; i < _rows; i++)
    {
        for (int j = 0; j < _cols; j++)
        {
            new_matrix._matrix[i][j] = _matrix[i][j] * other_matrix._matrix[i][j];
        }
    }
    return new_matrix;
}

template <typename T>
void Matrix_dense<T>::fill(const T &value)
{
    for (auto &row : _matrix)
    {
        std::fill(row.begin(), row.end(), value);
    }
}

template <typename T>
Matrix_dense<T> &Matrix_dense<T>::operator=(const T &value)
{
    fill(value);
    return *this;
}

template <typename T>
Matrix_dense<T> &Matrix_dense<T>::operator=(const Matrix_dense<T> &other_matrix)
{
    _rows = other_matrix._rows;
    _cols = other_matrix._cols;
    _matrix = other_matrix._matrix;

    return *this;
}

template <typename T>
bool Matrix_dense<T>::operator==(const Matrix_dense<T> &other_matrix) const
{
    if (_rows != other_matrix._rows || _cols != other_matrix._cols)
    {
        return false; // different dimensions
    }

    for (int i = 0; i < _rows; ++i)
    {
        if (!std::equal(_matrix[i].begin(), _matrix[i].end(), other_matrix._matrix[i].begin()))
        {
            return false;
        }
    }

    return true;
}

template <typename T>
Matrix_dense<T> Matrix_dense<T>::get_transpose()
{
    Matrix_dense<T> matrix_T(_cols, _rows);
    for (int j = 0; j < _cols; j++)
    {
        for (int i = 0; i < _rows; i++)
        {
            matrix_T._matrix[j][i] = _matrix[i][j];
        }
    }
    return matrix_T;
}

template <typename T>
void Matrix_dense<T>::transpose()
{
    *this = get_transpose();
}

template <typename T>
ostream &operator<<(ostream &output, const Matrix_dense<T> &matrix)
{
    for (int i = 0; i < matrix._rows; i++)
    {
        for (int j = 0; j < matrix._cols; ++j)
        {
            output << matrix._matrix[i][j] << " ";
        }
        output << endl;
    }
    return output;
}

template <typename T>
array<int, 2> Matrix_dense<T>::size() const
{
    return {_rows, _cols};
}

template <typename T>
int Matrix_dense<T>::nrows() const
{
    return _rows;
}

template <typename T>
int Matrix_dense<T>::ncols() const
{
    return _cols;
}

template <typename T>
T *Matrix_dense<T>::data()
{
    if (_contiguous_data.empty())
    {
        for (const auto &row : _matrix)
        {
            _contiguous_data.insert(_contiguous_data.end(), row.begin(), row.end());
        }
    }
    return _contiguous_data.data();
}

template <typename T>
const T *Matrix_dense<T>::data() const
{
    if (_contiguous_data.empty())
    {
        for (const auto &row : _matrix)
        {
            _contiguous_data.insert(_contiguous_data.end(), row.begin(), row.end());
        }
    }
    return _contiguous_data.data();
}

template <typename T>
void Matrix_dense<T>::resize(int new_rows, int new_cols)
{

    Matrix_dense<T> new_matrix(new_rows, new_cols) if (new_rows < _matrix._rows && new_cols < _matrix._cols)
    {
        for (int i = 0; i < new_rows; i++)
        {
            for (int j = 0; j < new_cols; j++)
            {
                new_matrix(i, j) = _matrix(i, j)
            }
        }
    }

    _matrix = new_matrix;
    _rows = new_rows;
    _cols = new_cols;
}

// template <typename T>
// T *Matrix_dense<T>::data()
// {
//     // if (_rows == 0 || _cols == 0)
//     return nullptr;
// }

// #include "lapack.h"

// template <typename T>
// std::tuple<Matrix_dense<T>, Matrix_dense<T>> solve_eigensystem(Matrix_dense<T> &A)
// {
//     // First, check if the matrix is symmetric
//     if (A.nrows() != A.ncols())
//     {
//         throw std::logic_error("solve_eigensystem only works with square matrices.");
//     }

//     // LAPACK expects double or float, make sure your Matrix_dense is specialized with those.
//     // Prepare the data for LAPACK++
//     int64_t n = A.nrows();
//     int64_t lda = n;
//     Matrix_dense<T> w(n, 1);        // Eigenvalues
//     Matrix_dense<T> work(1, 5 * n); // Workspace
//     int64_t lwork = 5 * n;
//     int64_t info;

//     lapack::syev('V', 'U', n, A.data(), lda, w.data(), work.data(), lwork, &info);

//     if (info != 0)
//     {
//         throw std::runtime_error("LAPACK execution failed.");
//     }

//     // Eigenvectors are now in A, eigenvalues in w.
//     return std::make_tuple(A, w);
// }

// template <typename T>
// void fill_random(std::vector<T> &vec, T min, T max)
// {
//     std::random_device rd;                     // Obtain a random number from hardware
//     std::mt19937 eng(rd());                    // Seed the generator
//     std::uniform_real_distribution<T> distr(); // Define the range

//     std::generate(vec.begin(), vec.end(), [&]()
//                   { return distr(eng); });
// }

// template <typename T>
// T lowest_eigenvalue(Matrix_sparse A, double tol)
// {
//     int n = A.nrows();
//     vector<T, n> v;
//     fill_random(v, 0, 1);

//     T normalizer = dot(v, v);
//     v = v.divide_by(normalizer);

//     vector<T> w = A.left_multiply(v);
//     T alpha = dot(v, w);
//     vector<T> v_old = 0; // how would i do this with template ??
//     T beta = zero;       // template???

//     vector<T, n> alphas;
//     vector<T, n - 1> betas;
//     Matrix_dense<T> TT;

//     for (int i = 0; i < n; i++)
//     {
//         w -= alpha * v + beta * v_old;
//         T w_normalizer = dot(w, w);
//         beta = w_normalizer;

//         if (beta < tol)
//             break;

//         v_old = v;
//         v = w / beta;
//         w = A.left_multiply(v);
//         T alpha_new = dot(v, w);

//         TT(i, i) = alpha;
//         if (i + 1 < n)
//         {
//             TT[i, i + 1] = beta;
//             TT[i + 1, i] = beta;
//         }
//     }

//     tuple<Matrix_dense<T>, Matrix_dense<T>> solve_eigensystem(TT);
//     /// okay print out and see if you need to sort to get the smallest or whats the deal
// }

// template <typename T>
// T davidson_method(Matrix_sparse &A)
// {

//     int n = A.nrows();
//     Matrix_dense<T> B(n, 0);
//     // add ability to add by column/ rows
//     vector<T> diags = A.get_diag();

//     vector<T, n> v;
//     fill_random(v, 0, 1);

//     T normalizer = dot(v, v);
//     v = v.divide_by(normalizer);

//     vector<T> v_old = v;
//     B.append(v);

//     vector<T> w;

//     for (int i = 0; i < n; i++)
//     {
//         w = A.left_multiply(v_old);

//         for (int j = 0; j < cols_; j++)
//         {
//         }
//     }
// }

template class Matrix_dense<double>;
