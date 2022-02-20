#pragma once

#include <vector>
#include <functional>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "Eigen/Dense"
#define M_PI 3.14159265358979323846

using namespace Eigen;

using condFunc = std::function<double(double)>;
using matrix = std::vector<std::vector<double>>;

struct ToCVector {
    matrix v1Ar;
    matrix v2Ar;
    matrix uAr;
};

struct coordProp {
    double x;
    int Nx;
    VectorXd xLin;
};



class DE
{
protected:
    double x;
    double y;
    double t;
    int Nx;
    int Ny;
    int Nt;
    double dx;
    double dy;
    double dt;

    condFunc stCond;
    condFunc leftCond;
    condFunc rightCond;
    condFunc y1Cond;
    condFunc y2Cond;

    VectorXd xLin;
    VectorXd tLin;
    VectorXd yLin;

    MatrixXd matU;

public:
    DE(double x,
        condFunc leftCond, condFunc rightCond,
        int Nx);
    virtual ~DE();

    MatrixXd GetU();
    coordProp getX();
    coordProp getY();
    coordProp getT();
    
    virtual MatrixXd startMatrix();

    static VectorXd ThomasAlg(MatrixXd& M, VectorXd& V);
    static ToCVector ar_cast(MatrixXd& m, VectorXd& Lin1, VectorXd& tLin2);
    friend std::ostream& operator<< (std::ostream& out, matrix& m);
    //static void printVectorMatrix(matrix m);
};

class WaveDE : public DE
{
private:
    double c;
    condFunc stCondDer;
public:
    WaveDE(double x, double t, double c,
        condFunc stCond, condFunc stCondDer,
        condFunc leftCond, condFunc rightCond,
        int Nx = 6, int Nt = 24);

    MatrixXd startMatrix() override;
    MatrixXd Explicit();
    MatrixXd Implicit();

};

class LaplasDE : public DE
{
public:
    LaplasDE(double x, double y,
        condFunc leftCond, condFunc rightCond,
        condFunc y1Cond, condFunc y2Cond,
        int Nx, int Ny);

    MatrixXd startMatrix() override;
    MatrixXd Solution(double eps, double omega, int& k);
    matrix k_from_omega();
    double calcChart(int i, int j);
};

class BEelement
{
private:
    double xl;
    double xr;
    unsigned int Nx;
    Matrix4d aMat;
    Matrix4d xMat;
    VectorXd xLin;

    RowVector4d PolyVector(double x);
    Matrix4d XMatrix();
    Matrix4d AMatrix();
public:
    BEelement(double x1, double x2);
    
    RowVector4d NVector(double x);
    std::vector<double> getX();
    std::vector<std::vector<double>> forPlot();
};