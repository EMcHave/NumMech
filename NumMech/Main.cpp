#include <string>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <memory>
#include "Eigen/Dense"
#include "FEM.h"
#include "NumMech.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;


int main()
{
    float l = 2.4;
    int n = 24;
    vector<Constraint> cons{ Constraint{0, 1, 0}, Constraint{n, 1, 0} };
    vector<Force> forces{ Force{ 6, -29419, 0 } };
    BeamFEM* beam = new BeamFEM(l, n, cons, true, true, pair<int, int>{6, 11}, forces);
    VectorXd disps = beam->Solve();
    cout << disps << endl;

    vector<float> X;
    vector<float> Y;
    for (int i = 0; i <= n; i++)
        X.push_back(beam->nodes[i].x);
    for(int i = 0; i < disps.size(); i+=2 )
        Y.push_back(disps[i]);

    plt::plot(X, Y);
    plt::title("Bended beam");
    plt::show();
    delete beam;
    return 0;
}

/*
int main(int argc, char** argv)
{
    TrussFEM* truss = new TrussFEM(argv[1]);
    VectorXd displacements = truss->Solve();
    VectorXd defs = truss->Deformations();
    VectorXd stresses = truss->Stresses();
    VectorXd forces = truss->Forces();

    cout << displacements << endl;

    ofstream VTK("Truss.vtk");
    ofstream difF("difF.txt");
    ofstream difU("difU.txt");
    ifstream abaqusF("forces.txt");
    ifstream abaqusU("disps.txt");
    ofstream U("dispTable.txt");
    ofstream F("forceTable.txt");

    VTK << "# vtk DataFile Version 1.0\n3D triangulation data\n" 
        << "ASCII\n\nDATASET POLYDATA\nPOINTS " << truss->nodes.size() << " float" << endl;

    for (const Node& n : truss->nodes)
    {
        VTK << n.x << ' ' << n.y << ' ' << 0 << endl;
    }

    VTK << "LINES " << truss->elements.size() << ' ' << 3 * truss->elements.size() << endl;

    for (Element* n : truss->elements)
    {
        TrussElement* N = dynamic_cast<TrussElement*>(n);
        VTK << 2 << '\t' << N->node_i.id << '\t' << N->node_j.id << endl;
    }
    
    VTK << "POINT_DATA " << truss->nodes.size() << endl;
    VTK << "VECTORS Displacements float" << endl;

    for (size_t i = 0; i < displacements.size(); i += 2)
    {
        VTK << displacements(i) << ' ' << displacements(i + 1) << ' ' << 0 << endl;
        U << i + 1  << '\t' << displacements(i) << '\t' << displacements(i + 1) << endl;
    }
    VTK << "CELL_DATA " << truss->elements.size() << endl;

    VTK << "SCALARS Forces float 1\nLOOKUP_TABLE default" << endl;
    for (size_t i = 0; i < forces.size(); i++)
    {
        VTK << fixed << forces(i) << endl;
        F << i+1 << '\t' << forces(i) << endl;
    }

    VectorXd abaqU(displacements.size()), abaqF(forces.size());
    int n,m, k;
    k = 0;
    while (!abaqusF.eof())
    {
        abaqusF >> n >> m >> abaqF(k++);
    }
    k = 0;
    while (!abaqusU.eof())
    {
        abaqusU >> n >> abaqU(k++) >> abaqU(k++);
    }
    
    for (size_t i = 0; i < forces.size(); i++)
        difF << 'F' << i + 1 << '\t' << abaqF(i) - forces(i) << endl;
    difF << "Difference in forces\n" << abaqF - forces << endl;
    difU << "Difference in displacements\n" << abaqU - displacements << endl;

    delete truss;
    
    return 0;
}
*/

/* 4 lab
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


/* 3 lab
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


/* 2 lab
int main()
{
    using namespace std;
    function stCond{
        [](double x) { return (x + 0.2) * sin(M_PI * x / 2); }
    };

    function stCondDer{
        [](double x) { return 1 + x * x; }
    };

    function leftCond{
        [](double t) {return  0; }
    };

    function rightCond{
        [](double t) {return 1.2*(t+1); }
    };

    WaveDE* de = new WaveDE(1, 0.5, 1, stCond, stCondDer, leftCond, rightCond, 6, 24);

    VectorXd xLin = de->getX().xLin;
    VectorXd tLin = de->getT().xLin;

    MatrixXd resultIm = de->Implicit();
    MatrixXd resultEx = de->Explicit();

    MatrixXd dif = resultEx - resultIm;

    ToCVector forPlot = DE::ar_cast(resultIm, xLin, tLin);

    ofstream file("outTable.txt");

    file << dif << endl;


    plt::plot_surface(forPlot.v1Ar, forPlot.v2Ar, forPlot.uAr, std::map<string, string>{ {"cmap", "plasma"}});
    plt::title("Implicit");
    plt::xlabel("x");
    plt::ylabel("t");
    plt::set_zlabel("U");
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
    

    delete de;

    return 0;
}
*/