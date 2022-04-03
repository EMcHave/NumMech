#include <string>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <memory>
#include "Eigen/Dense"
#include "FEM.h"
#include "BeamFEM.h"
#include "NumMech.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;

vector<MatrixXd> noAngles(MatrixXd&, MatrixXd&, MatrixXd&);


int main(int argc, int argv)
{
    float l = 2;
    int n = 20;
    vector<Constraint> cons{ Constraint{0, 1, 1}, Constraint{n, 1, 0} };
    vector<Force> forces{ Force{ 10, -98000, 0 } };
    BeamFEM* beam = new BeamFEM(l, n, cons, false, false, pair<int, int>{5, 9}, forces);

    vector<vector<VectorXd>> dynSol= beam->DynamicSolve();

    vector<float> X, Y1, Y2, Y3;
    VectorXd Xvec(n + 1);

    for (int i = 0; i <= n; i++)
    {
        X.push_back(beam->nodes[i].x/2);
        Xvec(i) = beam->nodes[i].x/2;
    }
        
    for (int i = 0; i < dynSol.at(0).at(0).size(); i += 2)
    {
        Y1.push_back(dynSol.at(0).at(100000)[i]); //100, 549.6, 99, 550, 99, 550
        Y2.push_back(dynSol.at(0).at(250000)[i]); //250, 700.1, 249, 700, 250.5, 700 
        Y3.push_back(dynSol.at(0).at(500000)[i]); //500, 989.8, 499, 990, 499.99, 989,5
    }

    int k = n + 1;

    MatrixXd U(2 * n + 2, 6), V(2 * n + 2, 6), A(2 * n + 2, 6);
    MatrixXd nU((n + 1), 6), nV((n + 1), 6), nA((n + 1), 6);

    ofstream displacements("dyndispC++.txt"), velocities("velocitiesC++.txt"), accelerations("accelC++.txt"), Xtable("X.txt");

    //ofstream dispsComp("dispsComp.txt"), anglesComp("anglesComp.txt"), forcesComp("beamforcesComp.txt"), momentsComp("beammomentsComp.txt");

    ifstream dispsAb("beamdispsBE.txt"), forcesAb("forcesbeam.txt");

    VectorXd dispV(k), anglesV(k), forcesV(k), momentsV(k), dispsCompV(k), forcesCompV(k), momentsCompV(k), dispsAbV(k), anglesAbV(k), forcesAbV(k), momentsAbV(k);

    VectorXd numbers = VectorXd::LinSpaced(n + 1, 1, n + 1);


    U << dynSol.at(0).at(100000), dynSol.at(0).at(250000), dynSol.at(0).at(500000), dynSol.at(0).at(550000), dynSol.at(0).at(700000), dynSol.at(0).at(990000);
    V << dynSol.at(1).at(101000), dynSol.at(1).at(251000), dynSol.at(1).at(499800), dynSol.at(1).at(550000), dynSol.at(1).at(699000), dynSol.at(1).at(991000);
    A << dynSol.at(2).at(100000), dynSol.at(2).at(250000), dynSol.at(2).at(500000), dynSol.at(2).at(550000), dynSol.at(2).at(700000), dynSol.at(2).at(989500);

    for (int i = 0; i < k; i++)
    {
        nU.row(i) = U.row(2 * i);
        nV.row(i) = V.row(2 * i);
        nA.row(i) = A.row(2 * i);
    }

    displacements << numbers << '\n' << nU;
    velocities << numbers << '\n' << nV;
    accelerations << numbers <<'\n'<< nA;
    Xtable << Xvec;

    
    plt::figure(0);
    plt::named_plot("t1 = 0,1 c", X, Y1);
    plt::named_plot("t2 = 0.25 c", X, Y2);
    plt::named_plot("t3 = 0.5 c", X, Y3);
    plt::grid(true);
    plt::xlabel("X, m");
    plt::ylabel("A, m/(s*s)");
    plt::title("A");
    plt::legend();
    plt::show();
    
    delete beam;
}


/*
int main()
{
    float l = 1;
    int n = 12;
    vector<Constraint> cons{ Constraint{0, 1, 0}, Constraint{n, 1, 0} };
    vector<Force> f{ Force{ 6, 0, 0 } };
    BeamFEM* beam = new BeamFEM(l, n, cons, true, true, pair<int, int>{3, 5}, f);
    
    VectorXd disps = beam->StaticSolve();
    cout << disps << endl;
    vector<pair<vector<double>, vector<double>>> plot = beam->ContinuousFunctions();

    vector<double> X;
    vector<double> Y;
    vector<double> Y1;
    VectorXd Xvec(n + 1);
    for (int i = 0; i <= n; i++)
    {
        X.push_back(beam->nodes[i].x);
        Xvec(i) = beam->nodes[i].x;
    }
        
    for(int i = 0; i < disps.size(); i += 2 )
        Y.push_back(disps[i]);
    for (int i = 1; i < disps.size(); i += 2)
        Y1.push_back(disps[i]);

    
    int k = beam->nodes.size();

    ofstream displacements("dispsC++.txt"), angles("anglesC++.txt"), forces("beamforcesC++.txt"), moments("beammomentsC++.txt"), Xtable("X.txt");

    ofstream dispsComp("dispsComp.txt"), anglesComp("anglesComp.txt"), forcesComp("beamforcesComp.txt"), momentsComp("beammomentsComp.txt");

    ifstream dispsAb("beamdispsBE.txt"), forcesAb("beamforces.txt");

    VectorXd dispV(k), anglesV(k), forcesV(k), momentsV(k), dispsCompV(k), forcesCompV(k), momentsCompV(k), dispsAbV(k), anglesAbV(k), forcesAbV(k), momentsAbV(k);

    VectorXd numbers = VectorXd::LinSpaced(n+1, 1, n+1);

    Xtable << Xvec;


    MatrixXd m(n+1, 3);
    MatrixXd m1(n + 1, 3);

    dispV = Eigen::Map<VectorXd, Eigen::Unaligned>(Y.data(), Y.size());
    anglesV = Eigen::Map<VectorXd, Eigen::Unaligned>(Y1.data(), Y1.size());

    m << numbers , dispV,  anglesV;
    displacements << m;



    forcesV = Eigen::Map< VectorXd, Eigen::Unaligned>(plot.at(2).second.data(), plot.at(2).second.size());
    momentsV = Eigen::Map< VectorXd, Eigen::Unaligned>(plot.at(1).second.data(), plot.at(1).second.size());

    m1 << numbers, forcesV, momentsV;

    forces << m1;
    
    float c, ch;
    int j = 0, t = 0;
    while (!dispsAb.eof())
        dispsAb >> c >> anglesAbV(j++) >> dispsAbV(t++);
    j = 0; t = 0;
    while (!forcesAb.eof())
        forcesAb >> c >> ch >> forcesAbV(j++) >> momentsAbV(t++);
    
    dispsComp << dispV - dispsAbV;
    anglesComp << anglesV - anglesAbV;
    forcesComp << -forcesV + forcesAbV;
    momentsComp << -momentsV + momentsAbV;
    

    plt::figure(0);
    plt::plot(X, Y);
    plt::xlabel("X, m");
    plt::ylabel("displacement, m");
    plt::title("Beam");
    plt::grid(true);

    plt::figure(1);
    plt::plot(plot[0].first, plot[0].second);
    plt::xlabel("X, m");
    plt::ylabel("Displacement, m");
    plt::title("Smooth beam");
    plt::grid(true);

    plt::figure(2);
    plt::plot(plot[1].first, plot[1].second);
    plt::xlabel("X, m");
    plt::ylabel("Moment, Nm");
    plt::title("Bending moments");
    plt::grid(true);

    plt::figure(3);
    plt::plot(plot[2].first, plot[2].second);
    plt::xlabel("X, m");
    plt::ylabel("Force, N");
    plt::title("Force");
    plt::grid(true);
    plt::show();
    
    delete beam;

    return 0;
}
*/


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