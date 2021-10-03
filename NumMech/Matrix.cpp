#pragma once

#include "Matrix.h"

using namespace std;

template <typename T>
Matrix<T>::Matrix(int i, int j, T& init) :r(i), c(j)
{
    mat = vector<vector<T>>(i, vector<T>(j, init));
}

template <typename T>
Matrix<T>::Matrix(int i, int j) :r(i), c(j)
{
    mat = vector<vector<T>>(i, vector<T>(j, 0));
}

template <typename T>
Matrix<T>::Matrix(vector<vector<T>>& m)
{
    for (int i = 0; i < r; i++)
        mat[i](m[i]);
    r = m.size();
    c = m[0].size();
}

template <typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& m)
{
    for (int i = 0; i < r; i++)
        mat.at(i)(m.get2DVector().at(i));

    return this;
}


template <typename T>
int Matrix<T>::cols() { return c; }

template <typename T>
int Matrix<T>::rows() { return r; }

template <typename T>
vector<T>& Matrix<T>::row(int i)
{
    return this->mat[i];
}

template <typename T>
vector<T> Matrix<T>::col(int j)
{
    vector<T> temp;
    for (int i = 0; i < r; i++)
        temp.push_back(mat[i][j]);
    return temp;
}

template <typename T>
vector<vector<T>> Matrix<T>::get2DVector()
{
    return mat;
}

template <typename T>
void Matrix<T>::print()
{
    for (vector vec : mat)
    {
        for (T el : vec)
            cout << std::setw(3) << el;
        cout << '\n';
    }
}

template <typename T>
T& Matrix<T>::operator() (int i, int j)
{
    return mat[i][j];
}


template <typename T>
Matrix<T>& Matrix<T>::operator*= (double mult)
{
    for (vector<T>& v : mat)
        for (T& num : v)
            num *= mult;
    return this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& m)
{
    Matrix<T> temp(r, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            temp(i, j) = mat[i][j] + m(i, j);
    return temp;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& m)
{
    Matrix<T> temp(r, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            temp(i, j) = mat[i][j] - m(i, j);
    return temp;
}

template <typename T>
Matrix<T>::~Matrix()
{ }




