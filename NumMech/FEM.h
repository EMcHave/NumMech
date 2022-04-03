#pragma once
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/SparseCholesky"


using namespace std;
using Eigen::VectorXd;
using Eigen::Triplet;
using Eigen::SparseMatrix;

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
public:
	int Id;
	void virtual CalculateStiffnessMatrix(vector<Triplet<double>>&) = 0;
	void virtual CalculateMassMatrix(vector<Triplet<double>>&) = 0;
	virtual ~Element();
};



class FEM
{
protected:
	vector<Constraint> constraints;
	vector<Force> forces;

	VectorXd forces_column;
	SparseMatrix<double> K_global;
	SparseMatrix<double> M_global;
	vector<Triplet<double>> triplets;
	vector<Triplet<double>> Mtriplets;

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


