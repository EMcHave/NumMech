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
private:
	const double E = 2*pow(10, 11);
	const double P = 0.3;
	const double A = 0.0001;

	array<int, 2> nodesIds;
public:
	Node node_i;
	Node node_j;
	int Id;
	double length;
	Element(Node, Node, int);
	void CalculateStiffnessMatrix(vector<Triplet<double>>&);
};

class TrussFEM
{
private:
	vector<Element> elements;
	vector<Node> nodes;
	vector<Constraint> constraints;
	vector<Force> forces;
	

	SparseMatrix<double> K_global;
	vector<Triplet<double>> triplets;

	VectorXd displacements;
	VectorXd defs;
	VectorXd stresses;
	VectorXd forces_in_truss;

	void readfile(const char*);
	void setConstraints(Eigen::SparseMatrix<double>::InnerIterator&, int);
	void applyConstraints();
public:
	TrussFEM(const char*);
	VectorXd Solve();
	VectorXd Deformations();
	VectorXd Stresses();
	VectorXd Forces();


};