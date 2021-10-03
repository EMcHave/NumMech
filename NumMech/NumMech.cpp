#pragma once

#include <iostream>
#include <algorithm>
#include <math.h>
#include "NumMech.h"


using condFunc = std::function<double(double)>;
using matrix = std::vector<std::vector<double>>;

using namespace Eigen;

/*****************************DE BASE CLASS*****************************/

DE::DE(double x, double t,
    condFunc stCond, 
    condFunc leftCond, condFunc rightCond, 
    int Nx, int Nt) 
    : x(x), t(t), Nt(Nt), Nx(Nx)
    {
        this->stCond = stCond;
        this->leftCond = leftCond;
        this->rightCond = rightCond;
        dx = x / Nx;
        dt = t / Nt;
        xLin = VectorXd::LinSpaced(Nx, 0, x);
        tLin = VectorXd::LinSpaced(Nt, 0, t);
    }

MatrixXd DE::GetU() { return matU; }

xProp DE::getX() {
    return xProp{ x, Nx, xLin };
}

tProp DE::getT() {
    return tProp{ t, Nt, tLin };
}

MatrixXd DE::startMatrix()
{
    MatrixXd matU = MatrixXd::Zero(Nt, Nx);
    matU.row(0) = xLin.unaryExpr(stCond);
    matU.col(0) = tLin.unaryExpr(leftCond);
    matU.col(Nx-1) = tLin.unaryExpr(rightCond);
    return matU;
}

VectorXd DE::ThomasAlg(MatrixXd &M, VectorXd &V)
{
    using namespace std;
    vector<double> delta, lambda;
    int s = V.size();

    delta.push_back(-M(0, 1) / M(0, 0));
    lambda.push_back(V(0) / M(0, 0));
    

    for (int i = 1; i < s-1; i++)
    {
        double denom = M(i, i - 1) * delta.at(i-1) + M(i, i);
        delta.push_back(-M(i, i + 1) / denom);
        lambda.push_back((V(i) - M(i, i - 1) * lambda.at(i-1)) / denom);
    }
    delta.push_back(0);
    lambda.push_back((V(s - 1) - M(s - 1, s - 2) * lambda.at(s - 2)) / (M(s - 1, s - 2) * delta.at(s-2) + M(s - 1, s - 1)));

    VectorXd resVector = VectorXd::Zero(s);
    resVector(s - 1) = lambda.at(s - 1);

    for (int i = s - 1; i > 0; i--)
        resVector(i - 1) = resVector(i) * delta.at(i - 1) + lambda.at(i - 1);

    return resVector;
}

ToCVector DE::ar_cast()
{
    matrix x(matU.rows());
    matrix t(matU.rows());
    matrix u(matU.rows());

    for (int i = 0; i < Nt; i++)
    {
        t[i] = std::vector<double>(Nx, tLin(i));
        x[i] = std::vector<double>(matU.cols());
        Map<RowVectorXd>(&x[i][0], 1, xLin.size()) = xLin;
        u[i] = std::vector<double>(matU.cols());
        Map<RowVectorXd>(&u[i][0], 1, matU.cols()) = matU.row(i);
    }

    return ToCVector{ x, t, u };
}

void DE::printVectorMatrix(matrix m)
{
    using namespace std;
    for (vector<double> vec : m)
    {
        for (double num : vec)
            cout << num << ' ';
        cout << '\n' << endl;
    }
}


/************************WAVEDE DERIVED CLASS**********************/


WaveDE::WaveDE(double x, double t, double c,
    condFunc stCond, condFunc stCondDer,
    condFunc leftCond, condFunc rightCond,
    int Nx, int Nt)
    :DE(x, t, stCond, leftCond, rightCond, Nx, Nt), c(c)
{
    this->dt = 0.1 * dx / c;
    this->stCondDer = stCondDer;
    matU = startMatrix();
}

MatrixXd WaveDE::startMatrix() 
{
    matU = DE::startMatrix();
    DE::startMatrix();
    VectorXd xLin = VectorXd::LinSpaced(Nx, 0, x);
    double coef = pow(c * dt / dx, 2);
    for (int i = 1; i < Nx - 1; i++)
        matU(1, i) = stCond(xLin(i)) + dt * stCondDer(xLin(i)) + 0.5*coef * (matU(0, i+1) - 2 * matU(0, i) + matU(0, i-1));

    return matU;
}

MatrixXd WaveDE::Explicit()
{
        
    double coef = pow(c * dt / dx, 2);

    for (int k = 1; k < Nt-1; k++)
        for (int i = 1; i < Nx-1; i++)
        {
            matU(k+1, i) = 2 * matU(k, i) - matU(k-1, i) + coef * (matU(k, i+1) - 2 * matU(k, i) + matU(k, i-1));
        }
 
    return matU;
}

MatrixXd WaveDE::Implicit()
{
    MatrixXd M = MatrixXd::Zero(Nx, Nx);
    M(0, 0) = -1;
    M(Nx - 1, Nx - 1) = -1;

    for (int i = 1; i < Nx-1; i++)
    {
        M(i, i - 1) = pow(dt, 2) / pow(dx, 2);
        M(i, i) = -(2 * pow(c * dt, 2) + pow(dx, 2)) / pow(dx * c, 2);
        M(i, i + 1) = M(i, i - 1);
    }

    for (int k = 1; k < Nt-1; k++) {
        VectorXd vec = 1 / pow(c, 2) * (-2 * matU.row(k) + matU.row(k - 1));
        vec(0) = -matU(k + 1, 0);
        vec(Nx-1) = -matU(k + 1, Nx-1);
        matU.row(k + 1) = DE::ThomasAlg(M, vec);
    }

    return matU;
}


