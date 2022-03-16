#include "BeamFEM.h"
#include "symbolicc++.h"
#include <chrono>

BeamElement::BeamElement(Node n1, Node n2, int id) : node_i(n1), node_j(n2)
{
    Id = id;
    nodesIds.at(0) = n1.id;
    nodesIds.at(1) = n2.id;
    length = abs(n1.x - n2.x);
    float l = length;
    K << 12 / (l * l), 6 / l, -12 / (l * l), 6 / l,
        6 / l, 4, -6 / l, 2,
        -12 / (l * l), -6 / l, 12 / (l * l), -6 / l,
        6 / l, 2, -6 / l, 4;
    M << 156, 22 * l, 54, -13 * l,
        22 * l, 4 * l * l, 13 * l, -3 * l * l,
        54, 13 * l, 156, -22 * l,
        -13 * l, -3 * l * l, -22 * l, 4 * l * l;
    K *= E * J / l;
    M *= ro * S * l / 420;

}

BeamElement::~BeamElement()
{
}


void BeamElement::CalculateStiffnessMatrix(vector<Triplet<double>>& triplets)
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            triplets.push_back(Triplet<double>(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 0, K(2 * i + 0, 2 * j + 0)));
            triplets.push_back(Triplet<double>(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 1, K(2 * i + 0, 2 * j + 1)));
            triplets.push_back(Triplet<double>(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 0, K(2 * i + 1, 2 * j + 0)));
            triplets.push_back(Triplet<double>(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 1, K(2 * i + 1, 2 * j + 1)));
        }
    }
}

void BeamElement::CalculateMassMatrix(vector<Triplet<double>>& triplets)
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            triplets.push_back(Triplet<double>(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 0, M(2 * i + 0, 2 * j + 0)));
            triplets.push_back(Triplet<double>(2 * nodesIds[i] + 0, 2 * nodesIds[j] + 1, M(2 * i + 0, 2 * j + 1)));
            triplets.push_back(Triplet<double>(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 0, M(2 * i + 1, 2 * j + 0)));
            triplets.push_back(Triplet<double>(2 * nodesIds[i] + 1, 2 * nodesIds[j] + 1, M(2 * i + 1, 2 * j + 1)));
        }
    }
}

void BeamFEM::readMesh(const char* path)
{
    string s;
    char c;
    ifstream inp(path);

    getline(inp, s);
    int i;
    while (s != "*Element")
    {
        Node newNode;
        inp >> i >> c >> newNode.x >> c >> newNode.y;
        inp.get();
        auto oldpos = inp.tellg();
        getline(inp, s);
        inp.seekg(oldpos);
        newNode.id = i - 1;
        nodes.push_back(newNode);
    }
    int n1, n2;
    getline(inp, s);
    while (s != "*End")
    {
        inp >> i >> c >> n1 >> c >> n2;
        inp.get();
        auto oldpos = inp.tellg();
        getline(inp, s);
        inp.seekg(oldpos);
        BeamElement* newEl = new BeamElement(nodes.at(n1 - 1), nodes.at(n2 - 1), i - 1);
        elements.push_back(newEl);
    }
}

BeamFEM::BeamFEM(float l, int n, vector<Constraint> cons, bool wf, bool qf, pair<int, int> sl, vector<Force> concF)
{
    qLength = sl;
    length = l;
    constraints = cons;
    forces = concF;
    
    float dx = l / n;
    for (int i = 0; i <= n; i++)
        nodes.push_back(Node{ i, i * dx, 0 });
    for (int i = 0; i < n; i++)
        elements.push_back(new BeamElement(nodes.at(i), nodes.at(i + 1), i));
    

    //readMesh("beamjob.inp");

    K_global = SparseMatrix<double>(2 * nodes.size(), 2 * nodes.size());
    M_global = SparseMatrix<double>(2 * nodes.size(), 2 * nodes.size());

    surface_forces = VectorXd::Zero(2 * nodes.size());
    vol_forces = VectorXd::Zero(2 * nodes.size());
    conc_forces = VectorXd::Zero(2 * nodes.size());




    for (Element* el : elements)
    {
        BeamElement* E = dynamic_cast<BeamElement*>(el);
        float dx = E->length;

        Eigen::Vector4f v;
        v << 1, dx / 6, 1, -dx / 6;
        Eigen::Vector4f vf = v * g * ro * S * dx / 2;
        Eigen::Vector4f sf = v * q * dx / 2;

        if (qf)
            if (E->Id >= sl.first && E->Id <= sl.second)
            {
                surface_forces(2 * E->node_i.id) += sf(0);
                surface_forces(2 * E->node_i.id + 1) += sf(1);
                surface_forces(2 * E->node_j.id) += sf(2);
                surface_forces(2 * E->node_j.id + 1) += sf(3);
            }

        if (wf)
        {
            vol_forces(2 * E->node_i.id) += vf(0);
            vol_forces(2 * E->node_i.id + 1) += vf(1);
            vol_forces(2 * E->node_j.id) += vf(2);
            vol_forces(2 * E->node_j.id + 1) += vf(3);
        }
    }
    for (vector<Force>::iterator it = forces.begin(); it != forces.end(); ++it)
    {
        conc_forces(2 * it->node_id) = it->fx;
        conc_forces(2 * it->node_id + 1) = it->fy;
    }

    forces_column = surface_forces + vol_forces + conc_forces;
}

VectorXd BeamFEM::StaticSolve()
{
    vector<Element*>::iterator it;

    for (it = elements.begin(); it != elements.end(); ++it)
    {
        (*it)->CalculateStiffnessMatrix(triplets);
    }

    K_global.setFromTriplets(triplets.begin(), triplets.end());

    applyConstraints();

    cout << Eigen::MatrixXd(K_global).determinant() << endl;


    //forces_column(2 * nodes.size() - 2) = 0.001;

    Eigen::BiCGSTAB<SparseMatrix<double>> solver;
    solver.setTolerance(pow(10, -15));
    solver.setMaxIterations(200);
    solver.compute(K_global);

    displacements = solver.solve(forces_column);

    return displacements;
}

vector<vector<VectorXd>> BeamFEM::DynamicSolve()
{

    int k = 1000000;
    float t = 1;
    float dt = t / k;

    vector<vector<VectorXd>> solution;
    vector<VectorXd> d2u(k), du(k), u(k);

    u.at(0) = Eigen::VectorXd::Zero(2 * nodes.size());
    u.at(1) = Eigen::VectorXd::Zero(2 * nodes.size());
    du.at(0) = Eigen::VectorXd::Zero(2 * nodes.size());
    d2u.at(0) = Eigen::VectorXd::Zero(2 * nodes.size());

    vector<Element*>::iterator it;

    for (it = elements.begin(); it != elements.end(); ++it)
    {
        (*it)->CalculateStiffnessMatrix(triplets);
        (*it)->CalculateMassMatrix(Mtriplets);
    }


    K_global.setFromTriplets(triplets.begin(), triplets.end());
    M_global.setFromTriplets(Mtriplets.begin(), Mtriplets.end());
    applyConstraints();

    Eigen::SimplicialLDLT<SparseMatrix<double>> explicit_solver;
    explicit_solver.compute(1 / (dt * dt) * M_global);
    

    Eigen::MatrixXd Kd = Eigen::MatrixXd(K_global);
    Eigen::MatrixXd Md = Eigen::MatrixXd(M_global);

    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 1; i < k - 1; i++)
    {
        VectorXd R = BeamFEM::CurrentForce(i * dt, t) * forces_column - (Kd - 2/(dt*dt)*Md) * u[i] - (1/(dt*dt)*Md)*u[i-1];
        u[i+1] = explicit_solver.solve(R);
        du[i] = 1 / (2 * dt) * (u[i + 1] - u[i - 1]);
        d2u[i] = 1 / (dt * dt) * (u[i + 1] - 2 * u[i] + u[i - 1]);
    }
    auto t2 = chrono::high_resolution_clock::now();
    auto ms_int = chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

    solution.push_back(u);
    solution.push_back(du);
    solution.push_back(d2u);

    return solution;
}

vector<pair<vector<double>, vector<double>>> BeamFEM::ContinuousFunctions()
{
    float l = length / elements.size();
    Symbolic a("a");
    Symbolic b("b");
    Symbolic k("k");
    Symbolic x("x");
    Symbolic Kappa, dKappa, V, N;
    VectorXd d = displacements;

    vector<double> v_field;
    vector<double> x_field;
    vector<double> x2_field;
    vector<double> M_field;
    vector<double> Q_field;

    N = (1 / (Symbolic(l * l * l)) * (x - b) * (2 * x * x - b * x - 3 * a * x - b * (b - 3 * a)),
        1 / (Symbolic(l * l)) * (x - a) * (x - b) * (x - b),
        1 / (Symbolic(l * l * l)) * (x - a) * (-2 * x * x + a * x + 3 * b * x + a * (a - 3 * b)),
        1 / (Symbolic(l * l)) * (x - b) * (x - a) * (x - a));

    Kappa = df(N, x, 2);
    dKappa = df(N, x, 3);


    for (int i = 0; i < elements.size(); i++)
    {
        BeamElement* el = dynamic_cast<BeamElement*>(elements[i]);
        double A = el->node_i.x;
        double B = el->node_j.x;
        double E = el->E;
        double J = el->J;
        V = (Symbolic(d(2 * i)), Symbolic(d(2 * i + 1)), Symbolic(d(2 * i + 2)), Symbolic(d(2 * i + 3)));
        Symbolic ABsubV = N[a == A, b == B];
        Symbolic ABsubM = Kappa[a == A, b == B];
        Symbolic ABsubQ = dKappa[a == A, b == B];
        Symbolic for_plotV = (V | ABsubV);
        Symbolic for_plotM = (V | ABsubM);
        Symbolic for_plotQ = (V | ABsubQ);
        VectorXd X = VectorXd::LinSpaced(5, A, B);
        for (int j = 0; j < X.size(); j++)
        {
            x_field.push_back(X(j));
            v_field.push_back(for_plotV[x == X(j)]);
            //M_field.push_back(-E*J*for_plotM[x == X(j)]);
            //Q_field.push_back(-E*J*for_plotQ[x == X(j)]);
        }

        x2_field.push_back(A);
        M_field.push_back(-E * J * for_plotM[x == A]);
        Q_field.push_back(-E * J * for_plotQ[x == A]);

        if (i == elements.size() - 1)
        {
            M_field.push_back(-E * J * for_plotM[x == B]);
            Q_field.push_back(-E * J * for_plotQ[x == B]);
            x2_field.push_back(B);
        }

    }

    vector<pair<vector<double>, vector<double>>> res;
    res.push_back(pair<vector<double>, vector<double>>(x_field, v_field));
    res.push_back(pair<vector<double>, vector<double>>(x2_field, M_field));
    res.push_back(pair<vector<double>, vector<double>>(x2_field, Q_field));
    return res;
}

BeamFEM::~BeamFEM()
{
    for (Element* e : elements)
        delete e;
}

float BeamFEM::CurrentForce(float t, float T)
{
    if (t <= T / 2)
        return 2 * t;
    else return 0;
}