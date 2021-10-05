#include <vector>
#include <string>
#include "Eigen/Dense"
#include "NumMech.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main()
{
    using namespace std;
    function stCond{
        [](double x) { return (x + 1) * sin(M_PI * x / 2); }
    };

    function stCondDer{
        [](double x) { return 1 - x*x; }
    };

    function leftCond{
        [](double t) {return  0.5*t; }
    };

    function rightCond{
        [](double t) {return 2; }
    };

    WaveDE* de = new WaveDE(1, 0.5, 1, stCond, stCondDer, leftCond, rightCond);

    VectorXd xLin = de->getX().xLin;
    VectorXd tLin = de->getT().tLin;

    MatrixXd resultIm = de->Implicit();
    MatrixXd resultEx = de->Explicit();

    MatrixXd dif = resultEx - resultIm;

    ToCVector forPlot = DE::ar_cast(dif, xLin, tLin);
    
    plt::plot_surface(forPlot.xAr, forPlot.tAr, forPlot.uAr);
    plt::xlabel("x");
    plt::ylabel("t");
    plt::set_zlabel("U");
    plt::show();


    
    delete de;
    
    return 0;
}
