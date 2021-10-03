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
        [](double x) { return (-x + 1) * cos(M_PI * x / 2); }
    };

    function stCondDer{
        [](double x) { return 2 * x + 1; }
    };

    function leftCond{
        [](double t) {return 2 * t + 1; }
    };

    function rightCond{
        [](double t) {return 0; }
    };

    WaveDE* de = new WaveDE(1, 0.5, 1, stCond, stCondDer, leftCond, rightCond);

    MatrixXd result = de->Implicit();

    ToCVector forPlot = de->ar_cast();
    
    plt::plot_surface(forPlot.xAr, forPlot.tAr, forPlot.uAr);
    plt::xlabel("x");
    plt::ylabel("t");
    plt::set_zlabel("U");
    plt::show();
    
    delete de;
    
    return 0;
}
