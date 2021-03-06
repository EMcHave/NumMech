#pragma once
#include "NumMech.h"


using condFunc = std::function<double(double)>;
using matrix = std::vector<std::vector<double>>;


/*****************************DE BASE CLASS*****************************/

DE::DE(double x, 
    condFunc leftCond, condFunc rightCond, 
    int Nx) 
    : x(x), Nx(Nx)
    {
        this->leftCond = leftCond;
        this->rightCond = rightCond;
        dx = x / Nx;
        xLin = VectorXd::LinSpaced(Nx, 0, x);
    }

DE::~DE(){}

MatrixXd DE::GetU() { return matU; }

coordProp DE::getX() {
    return coordProp{ x, Nx, xLin };
}

coordProp DE::getT() {
    return coordProp{ t, Nt, tLin };
}

coordProp DE::getY() {
    return coordProp{ y, Ny, yLin };
}

MatrixXd DE::startMatrix()
{
    MatrixXd temp = MatrixXd::Zero(Nt, Nx);
    temp.row(0) = xLin.unaryExpr(stCond);
    temp.col(0) = tLin.unaryExpr(leftCond);
    temp.col(Nx-1) = tLin.unaryExpr(rightCond);
    return temp;
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

ToCVector DE::ar_cast(MatrixXd &m, VectorXd &Lin1, VectorXd& Lin2)
{
    matrix v1(m.rows());
    matrix v2(m.rows());
    matrix u(m.rows());

    for (int i = 0; i < Lin2.size(); i++)
    {
        v2[i] = std::vector<double>(Lin1.size(), Lin2(i));
        v1[i] = std::vector<double>(m.cols());
        Eigen::Map<Eigen::RowVectorXd>(&v1[i][0], 1, Lin1.size()) = Lin2;
        u[i] = std::vector<double>(m.cols());
        Eigen::Map<Eigen::RowVectorXd>(&u[i][0], 1, m.cols()) = m.row(i);
    }

    return ToCVector{ v1, v2, u };
}

std::ostream& operator<<(std::ostream& out, matrix& m)
{
    using namespace std;
    for (vector<double> vec : m)
    {
        for (double num : vec)
            out << num << ' ';
        out << '\n' << endl;
    }
    return out;
}



/************************WAVEDE DERIVED CLASS**********************/


WaveDE::WaveDE(double x, double t, double c,
    condFunc stCond, condFunc stCondDer,
    condFunc leftCond, condFunc rightCond,
    int Nx, int Nt)
    :DE(x, leftCond, rightCond, Nx), c(c)
{
    this->Nt = Nt;
    this->t = t;
    this->stCond = stCond;
    this->stCondDer = stCondDer;

    dt = t / Nt;
    tLin = VectorXd::LinSpaced(Nt, 0, t);

}

MatrixXd WaveDE::startMatrix()
{
    MatrixXd temp = DE::startMatrix();

    double coef = pow(c * dt / dx, 2);
    for (int i = 1; i < Nx - 1; i++)
        temp(1, i) = stCond(xLin(i)) + dt * stCondDer(xLin(i)) + 0.5*coef * (temp(0, i+1) - 2 * temp(0, i) + temp(0, i-1));

    return temp;
}

MatrixXd WaveDE::Explicit()
{
    matU = WaveDE::startMatrix();
    double coef = (c * dt / dx) * (c * dt / dx);

    for (int k = 1; k < Nt-1; k++)
        for (int i = 1; i < Nx-1; i++)
        {
            matU(k+1, i) = 2 * matU(k, i) - matU(k-1, i) + coef * (matU(k, i+1) - 2 * matU(k, i) + matU(k, i-1));
        }
 
    return matU;
}

MatrixXd WaveDE::Implicit()
{;
    matU = WaveDE::startMatrix();
    
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



/*********************************** LAPLASDE class ***************************************/



LaplasDE::LaplasDE(double x, double y, condFunc leftCond, condFunc rightCond, condFunc y1Cond, condFunc y2Cond, int Nx, int Ny):
    DE(x, leftCond, rightCond, Nx)
{
    this->y = y;
    this->y1Cond = y1Cond;
    this->y2Cond = y2Cond;
    dy = dx;
    this->Ny = y/dy;
    yLin = VectorXd::LinSpaced(Ny, 0, y);

}

MatrixXd LaplasDE::Solution(double eps, double omega, int& k)
{
    matU = startMatrix();
    MatrixXd temp;
    do {
        temp = matU;

        for (int i = 1; i < Nx - 1; i++)
            for (int j = 1; j < Ny - 1; j++)
                matU(i, j) = matU(i, j) + omega * (calcChart(i, j) - matU(i, j));
        k++;
    } while ((matU - temp).lpNorm<Eigen::Infinity>() > eps);

    return matU;
}


matrix LaplasDE::k_from_omega()
{
    matrix temp(2);
    double omega = 0.1;
    int K = 0;
    while (omega < 2)
    {
        int k = 0;
        
        Solution(pow(10, -12), omega, k);
        temp[0].push_back(omega);
        temp[1].push_back(k);
        omega += 0.1;
        K += k;
    }
    std::cout << K << std::endl;
    return temp;
}


double LaplasDE::calcChart(int i, int j)
{
    return 0.25 * (matU(i - 1, j) + matU(i + 1, j) + matU(i, j - 1) + matU(i, j + 1));
}

MatrixXd LaplasDE::startMatrix()
{
    MatrixXd temp = MatrixXd::Zero(Nx, Ny);
    temp.row(0) = xLin.unaryExpr(y1Cond);
    temp.row(Ny - 1) = xLin.unaryExpr(y2Cond);
    temp.col(0) = yLin.unaryExpr(leftCond);
    temp.col(Nx - 1) = yLin.unaryExpr(rightCond);
    return temp;
}


//////////////////////// BEelement /////////////////////////


BEelement::BEelement(double x1, double x2): xl(x1), xr(x2), Nx(20)
{
    xMat = XMatrix();
    aMat = AMatrix();
    xLin = VectorXd::LinSpaced(Nx, xl, xr);
}

Matrix4d BEelement::XMatrix()
{
    Matrix4d temp;

    temp << 1, xl, xl* xl, xl* xl* xl,
        0, 1, 2 * xl, 3 * xl * xl,
        1, xr, xr* xr, xr* xr* xr,
        0, 1, 2 * xr, 3 * xr * xr;

    return temp;
}

Matrix4d BEelement::AMatrix()
{
    return xMat.inverse();
}

Eigen::RowVector4d BEelement::PolyVector(double x)
{
    return Eigen::RowVector4d{ 1, x, x * x, x * x * x };
}

Eigen::RowVector4d BEelement::NVector(double x)
{
    return PolyVector(x) * aMat;
}



std::vector<double> BEelement::getX()
{
    std::vector<double> temp;
    for (size_t i = 0; i < xLin.size(); i++)
        temp.push_back(xLin(i));
    return temp;
}

std::vector<std::vector<double>> BEelement::forPlot()
{
    std::ofstream out("funcForm.txt");

    std::vector<std::vector<double>> temp(4);
    for (size_t i = 0; i < xLin.size(); i++)
    {
        Eigen::RowVector4d N = NVector(xLin(i));
        out << std::setprecision(5) <<  N << '\n';
        temp.at(0).push_back(N(0));
        temp.at(1).push_back(N(1));
        temp.at(2).push_back(N(2));
        temp.at(3).push_back(N(3));
    }

    return temp;
}


