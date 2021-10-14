#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Dense"
#include "NumMech.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main()
{
    using namespace std;
    function y1Cond{
        [](double x) { return 10 * x * x * (1 - x); }
    };

    function y2Cond{
        [](double x) { return 50*sin(M_PI*x); }
    };

    function leftCond{
        [](double y) {return  30 * y * y * (1 - y); }
    };

    function rightCond{
        [](double y) {return 0; }
    };

    LaplasDE* de = new LaplasDE(1, 1, leftCond, rightCond, y1Cond, y2Cond, 20, 20);

    VectorXd xLin = de->getX().xLin;
    VectorXd yLin = de->getY().xLin;

    MatrixXd result = de->Solution(pow(10, -8), 1);


    ToCVector forPlot = DE::ar_cast(result, xLin, yLin);

    ofstream file("outTable.txt");

    
    
    plt::plot_surface(forPlot.v1Ar, forPlot.v2Ar, forPlot.uAr, std::map<string, string>{ {"cmap", "plasma"}});
    plt::title("Temperature");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("T");
    plt::show();
    
    /*
    plt::named_plot("First moment", forPlot.xAr[0], forPlot.uAr[0]);
    plt::named_plot("Second moment", forPlot.xAr[0], forPlot.uAr[15]);
    plt::named_plot("Third moment", forPlot.xAr[0], forPlot.uAr[45]);
    plt::title("Diffrent moments");
    plt::xlabel("x");
    plt::ylabel("U");
    plt::grid(true);
    plt::legend();
    //plt::set_zlabel("U");
    plt::show();
    */
    
    delete de;
    
    return 0;
}

