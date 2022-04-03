#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include "FEM.h"
using namespace std;



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
    
    for (int k = 0; k < M_global.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(M_global, k); it; ++it)
            for (vector<int>::iterator idit = indicesToConstraint.begin(); idit != indicesToConstraint.end(); ++idit)
                setConstraints(it, *idit);
                
}

Element::~Element()
{
}

FEM::~FEM()
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