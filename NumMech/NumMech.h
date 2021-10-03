#pragma once

#include <vector>
#include <functional>
#include "Eigen/Dense"
#define M_PI 3.14159265358979323846

using namespace Eigen;

using condFunc = std::function<double(double)>;
using matrix = std::vector<std::vector<double>>;

struct ToCVector {
    matrix xAr;
    matrix tAr;
    matrix uAr;
};

struct xProp {
    double x;
    int Nx;
    VectorXd xLin;
};

struct tProp {
    double t;
    int Nt;
    VectorXd tLin;
};
class DE
{
protected:
    double x;
    double t;
    int Nx;
    int Nt;
    double dx;
    double dt;
    condFunc stCond;
    condFunc leftCond;
    condFunc rightCond;
    MatrixXd matU;
    VectorXd xLin;
    VectorXd tLin;

public:
    DE(double x, double t,
        condFunc stCond,
        condFunc leftCond, condFunc rightCond,
        int Nx = 15, int Nt = 15);

    MatrixXd GetU();
    xProp getX();
    tProp getT();
    ToCVector ar_cast();

    virtual MatrixXd Explicit() = 0;
    virtual MatrixXd Implicit() = 0;
    virtual MatrixXd startMatrix();

    static VectorXd ThomasAlg(MatrixXd& M, VectorXd& V);
    static void printVectorMatrix(matrix m);
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
        int Nx = 25, int Nt = 25);

    MatrixXd startMatrix() override;
    MatrixXd Explicit() override;
    MatrixXd Implicit() override;

};

