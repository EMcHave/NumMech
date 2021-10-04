#pragma once
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <functional>

template <typename T>
class Matrix
{
private:
    std::vector<std::vector<T>> mat;
    //std::vector<std::vector<T>> matT;
    int c;
    int r;

public:
    Matrix(int i, int j, T& init);
    Matrix(int i, int j);
    Matrix(std::vector<std::vector<T>>& m);
    Matrix<T>& operator= (const Matrix<T>& m);

    int rows() const;
    int cols() const;
    std::vector<T>& row(int i);
    std::vector<T>& col(int j);
    std::vector<std::vector<T>> get2DVector() const;

    T& operator() (int i, int j);
    Matrix<T>& operator*= (double mult);
    Matrix<T> operator+ (const Matrix<T>& m);
    Matrix<T> operator- (const Matrix<T>& m);
    Matrix<T> Transpose();

    //void funcTo(std::function<double(double)> f);
    void funcToCol(std::function<double(double)> f, int j);
    void funcToRow(std::function<double(double)> f, int i);
    void transpose();
    void print();

    static std::vector<T> linspace(T start, T end, int num);
    template <typename T>
    friend std::ostream& operator<< (std::ostream& out, const Matrix<T>& m);
    template <typename T>
    friend std::ostream& operator<< (std::ostream& out, const std::vector<T>& vec);

    ~Matrix();
};


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
    for (int i = 0; i < m.size(); i++)
        mat.push_back(m[i]);
    r = m.size();
    c = m[0].size();
}

template <typename T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& m)
{
    for (int i = 0; i < r; i++)
        mat[i] = m.get2DVector()[i];

    return *this;
}



template <typename T>
int Matrix<T>::cols() const { return c; }

template <typename T>
int Matrix<T>::rows() const { return r; }

template <typename T>
vector<T>& Matrix<T>::row(int i)
{
    return this->mat[i];
}

template <typename T>
vector<T>& Matrix<T>::col(int j)
{
    *this = this->Transpose();
    return this->mat[j];
}

template <typename T>
vector<vector<T>> Matrix<T>::get2DVector() const
{
    return mat;
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
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& m)
{
    vector<vector<T>> temp(r, vector<T>(c));
    vector<vector<T>> temp2 = m.get2DVector();

    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            temp[i][j] = mat[i][j] + temp2[i][j];
    return Matrix<T>(temp);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& m)
{
    vector<vector<T>> temp(r, vector<T>(c));
    vector<vector<T>> temp2 = m.get2DVector();

    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            temp[i][j] = mat[i][j] - temp2[i][j];
    return Matrix<T>(temp);
}




template<typename T>
Matrix<T> Matrix<T>::Transpose()
{
    Matrix<T> temp(c, r);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            temp(j, i) = mat[i][j];
    return temp;
}

template<typename T>
void Matrix<T>::transpose()
{
    Matrix<T> temp = *this;
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            mat[j][i] = temp(i,j);
}




template<typename T>
void Matrix<T>::funcToCol(function<double(double)> f, int j)
{
    for (int t = 0; t < r; t++)
        mat[t][j] = f(mat[t][j]);

}

template<typename T>
void Matrix<T>::funcToRow(function<double(double)> f, int i)
{
    for (int t = 0; t < c; t++)
        mat[i][t] = f(mat[i][t]);
}

template <typename T>
static vector<T> Matrix<T>::linspace(T start, T end, int num)
{
    std::vector<double> linspaced;

    double start = static_cast<double>(start);
    double end = static_cast<double>(end);
    double num = static_cast<double>(num);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); 
    return linspaced;
}




template <typename T>
ostream& operator<< (ostream& out, const Matrix<T>& m)
{
    vector<vector<T>> mat = m.get2DVector();
 
    for (vector<T> vec : mat)
    {
        for (T el : vec)
            out << el << setw(10);
        out << '\n';
    }

    return out;
}

template <typename T>
ostream& operator<< (ostream& out, const vector<T>& vec)
{
    for (T el : vec)
        out << el << '\n';
    return out;
}

template <typename T>
void Matrix<T>::print()
{
    for (vector<T> vec : mat)
    {
        for (T el : vec)
            cout << el << std::setw(10);
        cout << '\n';
    }
}

template <typename T>
Matrix<T>::~Matrix()
{ }




