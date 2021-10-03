#include "Vector.h"

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