#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include "FEM.h"
using namespace std;
using namespace Eigen;

TrussElement::TrussElement(Node n1, Node n2, int id): node_i(n1), node_j(n2)
{
    Id = id;
	length = sqrt((n2.x - n1.x)* (n2.x - n1.x) + (n2.y - n1.y)* (n2.y - n1.y));
	nodesIds.at(0) = n1.id;
	nodesIds.at(1) = n2.id;
}

BeamElement::BeamElement(Node n1, Node n2, float l, int id) : node_i(n1), node_j(n2), length(l)
{
    Id = id;
    nodesIds.at(0) = n1.id;
    nodesIds.at(1) = n2.id;
    K << 12/(l*l), 6 / l, -12*(l*l), 6 / l,
        6 / l, 4, -6 / l, 2,
        -12/(l*l), -6 / l, 12/(l*l), -6 / l,
        6 / l, 2, -6 / l, 4;
    K *= E * J / l;
}


BeamElement::~BeamElement()
{
}

void TrussElement::CalculateStiffnessMatrix(vector<Triplet<double>>& triplets)
{
	Matrix2d B;
	Matrix4d K;
	Matrix<double, 2, 4> T;
    //cout << Id << ' ' << node_i.x << ' '<< node_i.y << ' ' << node_j.x << ' ' << node_j.y << ' ' << length << endl;
	B << 1, -1,
		-1,  1;
	B *= E * A / length;
	T << (node_j.x - node_i.x) / length, (node_j.y - node_i.y) / length, 0, 0,
		0, 0, (node_j.x - node_i.x) / length, (node_j.y - node_i.y)/length;
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


TrussElement::~TrussElement()
{
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
        (* it)->CalculateStiffnessMatrix(triplets);
    }

    K_global.setFromTriplets(triplets.begin(), triplets.end());
    
    applyConstraints();
    
    for (vector<Force>::iterator it = forces.begin(); it != forces.end(); ++it)
    {
        forces_column(2 * it->node_id) = it->fx;
        forces_column(2 * it->node_id + 1) = it->fy;
    }

    ConjugateGradient<SparseMatrix<double>> solver;
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
        double Y1 = itEl->node_i.y + displacements(2*itEl->node_i.id + 1);
        double X2 = itEl->node_j.x + displacements(2*itEl->node_j.id);
        double Y2 = itEl->node_j.y + displacements(2*itEl->node_j.id + 1);
        double newLength = sqrt((X2-X1)* (X2 - X1) + (Y2-Y1) * (Y2 - Y1));
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

TrussFEM::TrussFEM(const char* path)
{
    readfile(path);
    K_global = SparseMatrix<double>(2 * nodes.size(), 2 * nodes.size());
    forces_column = VectorXd::Zero(2 * nodes.size());
}

void FEM::setConstraints(SparseMatrix<double>::InnerIterator& it, int index)
{
    if (it.row() == index || it.col() == index)
    {
        it.valueRef() = it.row() == it.col() ? 1.0 : 0.0;
    }
}

void FEM::applyConstraints()
{
    vector<int> indicesToConstraint;

    for (vector<Constraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it)
    {
        if (it->x)
            indicesToConstraint.push_back(2 * it->node_id + 0);
        if (it->y)
            indicesToConstraint.push_back(2 * it->node_id + 1);
    }
    for (vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit)
        forces_column(*idit) = 0;

    for (int k = 0; k < K_global.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(K_global, k); it; ++it)
            for (vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit)
                setConstraints(it, *idit);
}

Element::~Element()
{
}

FEM::~FEM()
{
    for (auto* n : elements)
        delete n;
}

BeamFEM::BeamFEM(float l, int n, vector<Constraint> cons, bool wf, bool qf, pair<int, int> sl, vector<Force> concF)
{
    length = l;
    constraints = cons;
    forces = concF;
    float dx = l / n;
    for (int i = 0; i <= n; i++)
        nodes.push_back(Node{ i, i*dx, 0, 0 });
    for (int i = 0; i < n; i++)
        elements.push_back(new BeamElement(nodes.at(i), nodes.at(i + 1), dx, i));
    K_global = SparseMatrix<double>(2 * nodes.size(), 2 * nodes.size());
    surface_forces = VectorXd::Zero(2 * nodes.size());
    vol_forces = VectorXd::Zero(2 * nodes.size());
    conc_forces = VectorXd::Zero(2 * nodes.size());

;
    

    Vector4f v;
    v << 1, dx / 6, 1, -dx / 6;
    Vector4f vf = v * g * ro * S * dx / 2 * 0;
    Vector4f sf = v * q * dx / 2 ;

    for (Element* el : elements)
    {
        BeamElement* E = dynamic_cast<BeamElement*>(el);
        if (qf)
            if (E->Id >= sl.first && E->Id <= sl.second)
            {
                surface_forces(2 * E->node_i.id) = sf(0);
                surface_forces(2 * E->node_i.id + 1) = sf(1);
                surface_forces(2 * E->node_j.id) = sf(2);
                surface_forces(2 * E->node_j.id + 1) = sf(3);
            }

        if (wf)
        {
            vol_forces(2 * E->node_i.id) = vf(0);
            vol_forces(2 * E->node_i.id + 1) = vf(1);
            vol_forces(2 * E->node_j.id) = vf(2);
            vol_forces(2 * E->node_j.id + 1) = vf(3);
        }
    }
    for (vector<Force>::iterator it = forces.begin(); it != forces.end(); ++it)
    {
        conc_forces(2 * it->node_id) = it->fx;
        conc_forces(2 * it->node_id + 1) = it->fy;
    }

    forces_column = surface_forces + vol_forces + conc_forces;
}

VectorXd BeamFEM::Solve()
{
    vector<Element*>::iterator it;

    for (it = elements.begin(); it != elements.end(); ++it)
    {
        (*it)->CalculateStiffnessMatrix(triplets);
    }

    K_global.setFromTriplets(triplets.begin(), triplets.end());
    
    applyConstraints();

    ConjugateGradient<SparseMatrix<double>> solver;
    solver.setTolerance(pow(10, -15));
    solver.setMaxIterations(200);
    solver.compute(K_global);

    displacements = solver.solve(forces_column);

    return displacements;
}

BeamFEM::~BeamFEM()
{
}


/*
    set<int> loaded_nodes;
    for (Element* el : elements)
    {
        BeamElement* E = dynamic_cast<BeamElement*>(el);
        if (E->Id >= sl.first && E->Id <= sl.second)
            loaded_nodes.insert(E->node_i.id);
    }
    for (vector<Force>::iterator it = forces.begin(); it != forces.end(); ++it)
    {
        conc_forces(2 * it->node_id) = it->fx;
        conc_forces(2 * it->node_id + 1) = it->fy;
    }
    if (w)
    {

        for (int i = 0; i <= n; i+=2)
        {
            vol_forces(2 * nodes[i].id) = vf(0);
            vol_forces(2 * nodes[i].id + 1) = vf(1);
            vol_forces(2 * nodes[i+1].id + 0) = vf(2);
            vol_forces(2 * nodes[i+1].id + 1) = vf(3);
        }
    }
    //cout << vol_forces << '\n' << endl;


    for (set<int>::iterator it = loaded_nodes.begin(); it != loaded_nodes.end(); ++it)
    {

        surface_forces(4 * (*it)) = sf(0);
        surface_forces(4 * (*it) + 1) = sf(1);
        surface_forces(4 * (*it) + 2) = sf(2);
        surface_forces(4 * (*it) + 3) = sf(3);
    }
    forces_column = surface_forces + vol_forces + conc_forces;
*/