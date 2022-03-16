#pragma once
#include "FEM.h"


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
	void CalculateMassMatrix(vector<Triplet<double>>&) override;
	~TrussElement();
};


class TrussFEM : public FEM
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