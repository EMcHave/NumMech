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
	const double J = 3.5 * pow(10, -6);
	const double J1 = 6.42 * pow(10, -7);
	//const double J = 7.87*pow(10, -9);
	const double E = 2 * pow(10, 11);
	const float S = 0.000922;
	//const float S = 0.000144;
	const float ro = 7700;
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
	const float q = -3000;
	const float ro = 7700;
	const float S = 0.000922;
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
