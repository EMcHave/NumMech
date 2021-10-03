#pragma once

#include <vector>

template <typename T>
class Vector
{
public:
    Vector(int i);
    Vector(int i, T& init);
    Vector(std::vector<T>& vec);
    ~Vector();

private:
    std::vector<T> vec;
    int size;
};

template <typename T>
Vector<T>::Vector(int i) : vec(i) {}

template <typename T>
Vector<T>::Vector(int i, T& init) : vec(i, init) {}

template <typename T>
Vector<T>::Vector(std::vector<T>& vec)
{
    this->vec(vec);
}


template <typename T>
Vector<T>::~Vector()
{
    ~vec;
}