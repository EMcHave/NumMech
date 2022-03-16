#pragma once
#include "FEM.h"

class BeamElement : public Element
{
private:

	Eigen::Matrix4d K;
	Eigen::Matrix4d M;
	array<int, 2> nodesIds;
public:
	BeamElement(Node, Node, int);
	const double J = 4.0 * pow(10, -5);
	//double J = 5.09 * pow(10, -9);
	const double E = 2 * pow(10, 11);
	const float S = 0.0076;
	//float S = 0.000144;
	const float ro = 7900;
	float length;
	Node node_i;
	Node node_j;
	void CalculateStiffnessMatrix(vector<Triplet<double>>&) override;
	void CalculateMassMatrix(vector<Triplet<double>>&) override;
	~BeamElement();
};


class BeamFEM : public FEM
{
private:
	const float g = -9.81;
	const float q = -29420;
	const float ro = 7900;
	const float S = 0.0076;
	float length;
	VectorXd surface_forces;
	VectorXd vol_forces;
	VectorXd conc_forces;
	pair<int, int> qLength;

	void readMesh(const char* path);
public:
	BeamFEM(float l, int n, vector<Constraint>, bool wf, bool qf, pair<int, int>, vector<Force>);
	VectorXd StaticSolve();
	vector<vector<VectorXd>> DynamicSolve();
	vector<pair<vector<double>, vector<double>>> ContinuousFunctions();
	~BeamFEM();

	static float CurrentForce(float t, float T);
};
