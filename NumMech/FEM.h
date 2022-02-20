#pragma once
#include <vector>
#include <array>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"


using namespace std;
using namespace Eigen;

struct Node
{
	int id;
	double x;
	double y;
	double a;
};

struct Constraint
{
	int node_id;
	bool x;
	bool y;
};

struct Force
{
	int node_id;
	double fx;
	double fy;
};


class Element
{
public:
	int Id;
	void virtual CalculateStiffnessMatrix(vector<Triplet<double>>&) = 0;
	virtual ~Element();
};

class TrussElement : public Element
{
private: 
	const double E = 2 * pow(10, 11);
	const double P = 0.3;
	const double A = 0.0001;

	array<int, 2> nodesIds;
public:
	TrussElement(Node, Node, int);
	Node node_i;
	Node node_j;
	double length;
	void CalculateStiffnessMatrix(vector<Triplet<double>>&) override;
	~TrussElement();
};

class BeamElement : public Element
{
private:
	const double J = 5.027*pow(10, -9);
	const double E = 2 * pow(10, 11);
	Matrix4f K;
	array<int, 2> nodesIds;
public:
	BeamElement(Node, Node, float, int);
	float length;
	Node node_i;
	Node node_j;
	void CalculateStiffnessMatrix(vector<Triplet<double>>&) override;
	~BeamElement();
};


class FEM
{
protected:
	vector<Constraint> constraints;
	vector<Force> forces;

	VectorXd forces_column;
	SparseMatrix<double> K_global;
	vector<Triplet<double>> triplets;

	VectorXd displacements;
	VectorXd defs;
	VectorXd stresses;
	void setConstraints(Eigen::SparseMatrix<double>::InnerIterator&, int);
	void applyConstraints();
public: 
	vector<Element*> elements;
	vector<Node> nodes;
	virtual ~FEM();
};

class TrussFEM: public FEM
{
private:
	
	VectorXd forces_in_truss;

	void readfile(const char*);

public:

	TrussFEM(const char*);
	VectorXd Solve();
	VectorXd Deformations();
	VectorXd Stresses();
	VectorXd Forces();
	~TrussFEM();
};

class BeamFEM: public FEM
{
private:
	const float g = -9.8;
	const float q = -98066;
	const float ro = 7700;
	const float S = 1.6 * pow(10, -6);
	float length;
	VectorXd surface_forces;
	VectorXd vol_forces;
	VectorXd conc_forces;
public:
	BeamFEM(float l, int n, vector<Constraint>, bool wf, bool qf, pair<int, int>, vector<Force>);
	VectorXd Solve();
	VectorXd Deformations();
	VectorXd Stresses();
	VectorXd Forces();
	~BeamFEM();

};