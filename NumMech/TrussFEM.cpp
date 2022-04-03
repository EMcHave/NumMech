#include "TrussFEM.h"

TrussElement::TrussElement(Node n1, Node n2, int id) : node_i(n1), node_j(n2)
{
    Id = id;
    length = sqrt((n2.x - n1.x) * (n2.x - n1.x) + (n2.y - n1.y) * (n2.y - n1.y));
    nodesIds.at(0) = n1.id;
    nodesIds.at(1) = n2.id;
}

TrussElement::~TrussElement()
{
}


void TrussElement::CalculateStiffnessMatrix(vector<Triplet<double>>& triplets)
{
    Eigen::Matrix2d B;
    Eigen::Matrix4d K;
    Eigen::Matrix<double, 2, 4> T;
    //cout << Id << ' ' << node_i.x << ' '<< node_i.y << ' ' << node_j.x << ' ' << node_j.y << ' ' << length << endl;
    B << 1, -1,
        -1, 1;
    B *= E * A / length;
    T << (node_j.x - node_i.x) / length, (node_j.y - node_i.y) / length, 0, 0,
        0, 0, (node_j.x - node_i.x) / length, (node_j.y - node_i.y) / length;
    K = (T.transpose()) * B * T;
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

void TrussElement::CalculateMassMatrix(vector<Triplet<double>>&)
{
}


TrussFEM::TrussFEM(const char* path)
{
    readfile(path);
    K_global = SparseMatrix<double>(2 * nodes.size(), 2 * nodes.size());
    forces_column = VectorXd::Zero(2 * nodes.size());
}


void TrussFEM::readfile(const char* path)
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
    while (s != "*Constraints")
    {
        inp >> i >> c >> n1 >> c >> n2;
        inp.get();
        auto oldpos = inp.tellg();
        getline(inp, s);
        inp.seekg(oldpos);
        TrussElement* newEl = new TrussElement(nodes.at(n1 - 1), nodes.at(n2 - 1), i - 1);
        elements.push_back(newEl);
    }
    getline(inp, s);
    while (s != "*Forces")
    {
        Constraint newCon;
        inp >> i >> c >> newCon.x >> c >> newCon.y;
        newCon.node_id = i - 1;
        inp.get();
        auto oldpos = inp.tellg();
        getline(inp, s);
        inp.seekg(oldpos);
        constraints.push_back(newCon);
    }
    getline(inp, s);
    while (s != "*End")
    {
        Force newF;
        inp >> i >> c >> newF.fx >> c >> newF.fy;
        newF.node_id = i - 1;
        inp.get();
        auto oldpos = inp.tellg();
        getline(inp, s);
        inp.seekg(oldpos);
        forces.push_back(newF);
    }
}

VectorXd TrussFEM::Solve()
{
    vector<Element*>::iterator it;

    for (it = elements.begin(); it != elements.end(); ++it)
    {
        (*it)->CalculateStiffnessMatrix(triplets);
    }

    K_global.setFromTriplets(triplets.begin(), triplets.end());

    applyConstraints();

    for (vector<Force>::iterator it = forces.begin(); it != forces.end(); ++it)
    {
        forces_column(2 * it->node_id) = it->fx;
        forces_column(2 * it->node_id + 1) = it->fy;
    }

    Eigen::ConjugateGradient<SparseMatrix<double>> solver;
    solver.setTolerance(pow(10, -15));
    solver.setMaxIterations(200);
    solver.compute(K_global);

    displacements = solver.solve(forces_column);

    return displacements;
}

VectorXd TrussFEM::Deformations()
{
    defs = VectorXd(elements.size());
    for (vector<Element*>::iterator itE = elements.begin(); itE != elements.end(); ++itE)
    {
        TrussElement* itEl = dynamic_cast<TrussElement*>(*itE);
        double X1 = itEl->node_i.x + displacements(2 * itEl->node_i.id);
        double Y1 = itEl->node_i.y + displacements(2 * itEl->node_i.id + 1);
        double X2 = itEl->node_j.x + displacements(2 * itEl->node_j.id);
        double Y2 = itEl->node_j.y + displacements(2 * itEl->node_j.id + 1);
        double newLength = sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1));
        defs(itEl->Id) = (newLength - itEl->length) / itEl->length;
    }
    return defs;
}

VectorXd TrussFEM::Stresses()
{
    stresses = defs * 2 * pow(10, 11);
    return stresses;
}

VectorXd TrussFEM::Forces()
{
    forces_in_truss = stresses * 0.0001;
    return forces_in_truss;
}

TrussFEM::~TrussFEM()
{
}
