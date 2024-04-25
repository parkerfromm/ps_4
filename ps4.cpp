#include <vector>
#include <iostream>
#include <algorithm>
#include <random>
#include <array>
#include <tuple>
#include <memory>
#include <cmath> // for sqrt function
#include "myclasses.h"
#include "lapacke.h" // Correct include for LAPACK
using namespace std;

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

//     // lapack::syev('V', 'U', n, A.data(), lda, w.data(), work.data(), lwork, &info);
//     LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, A.data(), lda, w.data());

//     if (info != 0)
//     {
//         throw std::runtime_error("LAPACK execution failed.");
//     }

//     // Eigenvectors are now in A, eigenvalues in w.
//     return std::make_tuple(A, w);
// }

// template <typename T>
// T dot(vector<T> a, vector<T> b)
// {

//     T sum = 0;
//     for (int i = 0; i < a.size(); i++)
//     {
//         sum += a[i] * b[i]
//     }

//     return sum;
// };

// template <typename T>
// T dot(const std::vector<T> &a, const std::vector<T> &b)
// {
//     T sum = T(0);
//     for (size_t i = 0; i < a.size(); ++i)
//     {
//         sum += a[i] * b[i];
//     }
//     return sum;
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

template <typename T>
void fill_random(std::vector<T> &vec, T min, T max)
{
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<T> distr(min, max);

    for (auto &element : vec)
    {
        element = distr(eng);
    }
}

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

//     cout << "before lapack" << endl;
//         for (int i = 0; i < TT.nrows(); i++)
//         {
//             for (int j = 0; j < TT.ncols(); ++j)
//             {
//                 cout << TT(i, j) << " ";
//             }
//             cout << endl;
//         }
//         cout<<TT << endl;

//     tuple<Matrix_dense<T>, Matrix_dense<T>> solve_eigensystem(TT);
//     /// okay print out and see if you need to sort to get the smallest or whats the deal
// }

template <typename T>
T dot(const std::vector<T> &a, const std::vector<T> &b)
{
    T sum = T(0);
    for (size_t i = 0; i < a.size(); ++i)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

template <typename T>
T *convert_to_array(Matrix_dense<T> A)
{
    T *array = new T[A.nrows() * A.ncols()];
    for (size_t i = 0; i < A.nrows(); ++i)
    {
        for (size_t j = 0; j < A.ncols(); ++j)
        {
            array[i * A.ncols() + j] = A(i, j);
        }
    }
    return array;
}

template <typename T>
std::pair<std::vector<T>, Matrix_dense<T>> solve_eigensystem(Matrix_dense<T> A)
{
    int n = static_cast<int>(A.nrows());
    int lda = n;
    int info;

    // Allocate memory using smart pointers to manage life-cycle
    std::unique_ptr<T[]> a(convert_to_array(A)); // Matrix data
    std::unique_ptr<T[]> w(new T[n]);            // Eigenvalues

    // Solve eigenproblem
    info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, a.get(), lda, w.get());

    // Check for successful completion
    if (info > 0)
    {
        std::cerr << "The algorithm failed to compute eigenvalues. Error code: " << info << std::endl;
        throw std::runtime_error("Failed to compute eigenvalues.");
    }

    // Construct eigenvalue vector and eigenvector matrix
    std::vector<T> eigenvalues(w.get(), w.get() + n);
    Matrix_dense<T> eigenvectors(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            eigenvectors(i, j) = a[i * lda + j];
        }
    }

    return {eigenvalues, eigenvectors};
}

template <typename T>
std::pair<std::vector<T>, Matrix_dense<T>> lowest_eigenvalue(const Matrix_sparse<T> A, int max_it, double tol)
{
    int n = A.nrows();
    std::vector<T> v(n), w(n), v_old(n);
    // fill_random(v, T(0), T(1));
    std::default_random_engine generator;
    std::normal_distribution<T> distribution(-1.0, 2.0);
    std::vector<T> alphas, betas;

    for (size_t i = 0; i < n; ++i)
    {
        v[i] = distribution(generator);
    }

    T normalizer = sqrt(dot(v, v));
    for (auto &element : v)
        element /= normalizer;

    w = A.left_multiply(v);
    T alpha = dot(v, w);
    std::fill(v_old.begin(), v_old.end(), T(0));
    T beta = T(0);
    betas.push_back(beta);

    Matrix_dense<T> TT(max_it, max_it); // Assuming enough memory for a dense matrix
                                        // Initialize to zero

    int its = 0;
    for (int i = 0; i < max_it; ++i)
    {
        w = A.left_multiply(v);
        alpha = dot(v, w);
        alphas.push_back(alpha);
        for (int j = 0; j < n; j++)
            w[j] -= alphas[i] * v[j] + betas[i] * v_old[j];
        normalizer = sqrt(dot(w, w));
        beta = normalizer;

        betas.push_back(beta);

        if (betas[i + 1] < tol)
        {
            its = i + 1;
            std::cout
                << "Convergence reached at iteration " << its << std::endl;
            alphas.resize(i + 1);
            betas.resize(i + 2); // Include last valid beta
            TT.resize(i + 1, i + 1);
            TT.fill(T(0)); // Adjust size to actual computation
            break;
        } // Convergence check

        v_old = v;
        if (betas[i + 1] != 0)
        {
            for (int j = 0; j < n; j++)
                v[j] = w[j] / betas[i + 1];
        }
    }
    // v_old = v;
    // for (auto &element : v)
    //     element = w[i] / beta;
    // w = A.left_multiply(v);
    // alpha = dot(v, w);

    for (int i = 0; i < its; i++)
    {
        TT(i, i) = alphas[i];
        if (i > 0)
        {
            TT(i, i - 1) = betas[i];
            TT(i - 1, i) = betas[i];
        }
    }

    // TT(i, i) = alphas[i];
    // if (i > 0)
    // {
    //     TT(i, i - 1) = betas[i];
    //     TT(i - 1, i) = betas[i];
    // }

    // v_old = v;
    // if (betas[i + 1] != 0)
    // {
    //     for (int j = 0; j < n; j++)
    //         v[j] = w[j] / betas[i + 1];
    // }

    cout << "before lapack" << endl;
    for (int i = 0; i < TT.nrows(); i++)
    {
        for (int j = 0; j < TT.ncols(); ++j)
        {
            cout << TT(i, j) << " ";
        }
        cout << endl;
    }
    // cout<<TT << endl;

    //     auto [eigenvectors, eigenvalues] = solve_eigensystem(TT); // Assuming this returns sorted eigenvalues
    //     return eigenvalues(0, 0);                                 // Return the lowest eigenvalue

    std::vector<T> eigenvalues;
    Matrix_dense<T> eigenvectors;
    std::tie(eigenvalues, eigenvectors) = solve_eigensystem(TT);

    return {eigenvalues, eigenvectors};
};

// template <typename T>
// std::pair<std::vector<T>, Matrix_dense<T>> lowest_eigenvalue(const Matrix_sparse<T> A, int max_it, double tol)
// {
//     int n = A.nrows();
//     std::vector<T> v(n), w(n), v_old(n);
//     // fill_random(v, T(0), T(1));
//     std::default_random_engine generator;
//     std::normal_distribution<T> distribution(0.0, 1.0);
//     std::vector<T> alphas, betas;

//     for (size_t i = 0; i < n; ++i)
//     {
//         v[i] = distribution(generator);
//     }

//     T normalizer = sqrt(dot(v, v));
//     for (auto &element : v)
//         element /= normalizer;

//     w = A.left_multiply(v);
//     T alpha = dot(v, w);
//     std::fill(v_old.begin(), v_old.end(), T(0));
//     T beta = T(0);

//     Matrix_dense<T> TT(max_it, max_it); // Assuming enough memory for a dense matrix
//     TT.fill(T(0));

//     for (int i = 0; i < max_it; ++i)
//     {
//         w = A.left_multiply(v);
//         alpha = dot(v, w);
//         alphas.push_back(alpha);
//         for (int j = 0; j < n; j++)
//             w[j] -= alpha * v[j] + beta * v_old[j];
//         normalizer = sqrt(dot(w, w));
//         beta = normalizer;

//         betas.push_back(beta);

//         if (betas[i + 1] < tol)
//         {
//             std::cout
//                 << "Convergence reached at iteration " << i + 1 << std::endl;
//             // alpha.resize(j + 1);
//             // beta.resize(j + 2);        // Include last valid beta
//             TT.resize(i + 1, i + 1); // Adjust size to actual computation
//             break;
//         } // Convergence check

//         v_old = v;
//         for (auto &element : v)
//             element = w[i] / beta;
//         w = A.left_multiply(v);
//         alpha = dot(v, w); // Initialize to zero

//         TT(i, i) = alpha;
//         if (i > 0)
//         {
//             TT(i, i - 1) = beta;
//             TT(i - 1, i) = beta;
//         }

//         v_old = v;
//         if (beta != 0)
//         {
//             for (int j = 0; j < n; j++)

//                 v[j] = w[j] / beta;
//         }
//     }
//     cout << "before lapack" << endl;
//     for (int i = 0; i < TT.nrows(); i++)
//     {
//         for (int j = 0; j < TT.ncols(); ++j)
//         {
//             cout << TT(i, j) << " ";
//         }
//         cout << endl;
//     }
//     // cout<<TT << endl;

//     // auto [eigenvectors, eigenvalues] = solve_eigensystem(TT); // Assuming this returns sorted eigenvalues
//     // return eigenvalues(0, 0);

//     std::vector<T> eigenvalues;
//     Matrix_dense<T> eigenvectors;
//     std::tie(eigenvalues, eigenvectors) = solve_eigensystem(TT);

//     return {eigenvalues, eigenvectors};
// }

template <typename T>
void initialize_symmetric_tridiagonal(Matrix_sparse<T> &matrix, T b, T c)
{
    // matrix.resize(n, n);
    int n = matrix.nrows();
    for (int i = 0; i < n; i++)
    {
        matrix(i, i) = b;
        if (i > 0)
        {
            matrix(i, i - 1) = c;
            matrix(i - 1, i) = c;
        }
    }
};
//
// template <typename T>
// bool test_lowest_eigenvalue(int n)
// {
//     Matrix_sparse<T> A(n, n);
//     initialize_symmetric_tridiagonal(A, T(2), T(-1));

//     std::pair<std::vector<T>, Matrix_dense<T>> computed= lowest_eigenvalue(A);
//     cout<<"Eigenvalues: { "
//     for (int i = 0; i < computged.first.size(); i++){
//         cout << computed.first[i] << " ";
//     }
//     cout << endl;
//     T expected_lowest_eigenvalue = 2 - 2 * cos(M_PI / (n + 1)); // Analytical eigenvalue for this matrix
//     std::cout << "Computed Lowest Eigenvalue: " << computed_lowest_eigenvalue << std::endl;
//     std::cout << "Expected Lowest Eigenvalue: " << expected_lowest_eigenvalue << std::endl;
//     return std::abs(computed_lowest_eigenvalue - expected_lowest_eigenvalue) < 1e-5; // Tolerance
// }

template <typename T>
T davidson_method(Matrix_sparse<T> &A, int max_iterations, T tolerance)
{

    int n = A.nrows();
    Matrix_dense<T> B(n, n);
    B.fill(T(0));
    // add ability to add by column/ rows
    vector<T> diags = A.get_diag();

    vector<T> v(n);
    fill_random(v, T(0), T(1));

    T normalizer = sqrt(dot(v, v));
    for (auto &element : v)
        element /= normalizer;

    vector<T> v_old = v;
    B.append_column(v);

    vector<T> w(n);
    T alpha;

    for (int i = 0; i < max_iterations; i++)
    {
        w = A.left_multiply(v_old);
        alpha = dot(v_old, w) / dot(v_old, v_old);
        std::vector<T> residual = w;

        for (int j = 0; j < n; j++)
        {
            residual[j] -= alpha * v_old[j];
        }

        T res_norm = sqrt(dot(residual, residual));
        if (res_norm < tolerance)
        {
            std::cout << "Convergence reached at iteration " << (i + 1) << std::endl;
            break;
        }

        std::vector<T> correction(n);
        for (int j = 0; j < n; j++)
        {
            if (std::abs(diags[j] - alpha) > tolerance)
            {
                correction[j] = residual[j] / (diags[j] - alpha);
            }
            else
            {
                correction[j] = residual[j] / tolerance; // Avoid division by zero or very small numbers
            }
        }

        normalizer = sqrt(dot(correction, correction));
        for (auto &element : correction)
            element /= normalizer;

        v_old = correction;
        B.append_column(correction); // Update the subspace
    }

    return alpha;
};

int main()
{
    int n = 10; // Matrix size
    double tolerance = 1e-6;
    int max_iterations = 100;

    Matrix_sparse<double> A(n, n);
    initialize_symmetric_tridiagonal(A, 2.0, -1.0);
    double expected_lowest_eigenvalue = 2 - 2 * cos(M_PI / (n + 1));
    std::pair<std::vector<double>, Matrix_dense<double>> computed = lowest_eigenvalue(A, max_iterations, tolerance);

    cout << "Eigenvalues: { ";
    for (int i = 0; i < computed.first.size(); i++)
    {
        cout << computed.first[i] << " ";
    }
    cout << "}" << endl;

    std::vector<double> eigenvalues;
    Matrix_dense<double> eigenvectors;
    Matrix_dense<double> dense_A = A.convert_to_dense();
    cout << "before lapack" << endl;
    for (int i = 0; i < dense_A.nrows(); i++)
    {
        for (int j = 0; j < dense_A.ncols(); ++j)
        {
            cout << dense_A(i, j) << " ";
        }
        cout << endl;
    }

    std::tie(eigenvalues, eigenvectors) = solve_eigensystem(dense_A);

    cout << "expected Eigenvalues: { ";
    for (int i = 0; i < eigenvalues.size(); i++)
    {
        cout << eigenvalues[i] << " ";
    }
    cout << "}" << endl;

    double lowest_eig = davidson_method(A, 100, tolerance);
    cout << "lowest eig =" << lowest_eig << endl;

    cout << "expected lowest eig =" << expected_lowest_eigenvalue << endl;
}

// the code is not correct. I do not know if I am fillinf in
// my sparse matrix class correctly
// i should get a tridiagnol matrix but i am only getting a diagnol one when converting to dense.