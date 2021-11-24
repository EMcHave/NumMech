#include <string>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <memory>
#include "Eigen/Dense"
#include "FEM.h"
#include "NumMech.h"
//#include "matplotlibcpp.h"

//namespace plt = matplotlibcpp;
using namespace std;

int main()
{
    TrussFEM* truss = new TrussFEM("Truss.txt");
    VectorXd displacements = truss->Solve();
    VectorXd defs = truss->Deformations();
    VectorXd stresses = truss->Stresses();
    VectorXd forces = truss->Forces();

    cout << displacements << '\n' << endl;
    cout << defs << '\n' << endl;
    cout << stresses << '\n' << endl;
    cout << forces << '\n' << endl;
    delete truss;
    
    return 0;
}


/*
int main()
{
    double x1 = -1;
    double x2 = 1;
    
    unique_ptr<BEelement> elt(new BEelement(-1, 1));

    vector<vector<double>> forPlot = elt->forPlot();
    vector<double> xLin = elt->getX();


    vector<double> sum;
    for (size_t i = 0; i < forPlot.at(0).size(); i++)
        sum.push_back(forPlot.at(0).at(i) + forPlot.at(2).at(i));

    plt::plot(xLin, sum);
    plt::title("Third property of form funcs");
    plt::grid(true);
    plt::show();

    plt::named_plot("Displacemnt in left node", xLin, forPlot.at(0));
    plt::named_plot("Rotation in left node", xLin, forPlot.at(1));
    plt::named_plot("Displacemnt in right node", xLin, forPlot.at(2));
    plt::named_plot("Rotation in right node", xLin, forPlot.at(3));
    plt::xlabel("x");
    plt::ylabel("Value");
    plt::grid(true);
    plt::legend();
    plt::show();
}
*/


/*
int main()
{
    using namespace std;
    int k = 0;
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

    LaplasDE* de = new LaplasDE(1, 1, leftCond, rightCond, y1Cond, y2Cond, 10, 10);

    VectorXd xLin = de->getX().xLin;
    VectorXd yLin = de->getY().xLin;

    MatrixXd result = de->Solution(pow(10, -10), 1.5, k);

    auto t1 = std::chrono::high_resolution_clock::now();
    matrix k_omega = de->k_from_omega();
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    cout << duration << endl;


    ToCVector forPlot = DE::ar_cast(result, xLin, yLin);

    ofstream file("outTable.txt");

    for (vector vec : k_omega) {
        for (double d : vec)
            file << setw(5) << d;
        file << '\n' << endl;
    }
        


    
    plt::plot_surface(forPlot.v1Ar, forPlot.v2Ar, forPlot.uAr, std::map<string, string>{ {"cmap", "plasma"}});
    plt::title("Measured value");
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("U");
    plt::show();
    
    plt::named_plot("y1", forPlot.v1Ar[0], forPlot.uAr[0]);
    plt::named_plot("y2 > y1", forPlot.v1Ar[0], forPlot.uAr[5]);
    plt::named_plot("y3 > y2", forPlot.v1Ar[0], forPlot.uAr[9]);
    plt::title("Diffrent layers");
    plt::xlabel("x");
    plt::ylabel("U");
    plt::grid(true);
    plt::legend();
    plt::show();
    
    plt::plot(k_omega[0], k_omega[1]);
    plt::xlabel("Omega");
    plt::ylabel("Iterations");
    plt::grid(1);

    plt::show();
    
   
    
    delete de;
    
    return 0;
}
*/