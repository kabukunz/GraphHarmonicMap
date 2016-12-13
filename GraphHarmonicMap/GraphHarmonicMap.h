#pragma once

#include "Mesh/Mesh.h"
#include "Mesh/Iterators.h"
#include "Graph.h"

using namespace MeshLib;

#ifndef EPS
#define EPS 1e-8
#endif // !EPS


/*
* represent a target on the graph, identified by an arc, and the direction, and the length from starting point
* for a loop, there is only one starting point, but two directions: +1, -1;
* for otherwise, there are two starting points +1 and -1, direction +1 means starting from +1, -1 starting from -1
*/
class CTarget
{
public:
    SmartGraph::Edge edge;
    SmartGraph::Node node;
    double length;
};

typedef map<int, vector<int>> Cut;

class CGraphHarmonicMap
{
public:
    CGraphHarmonicMap();
    ~CGraphHarmonicMap();

    int setMesh(const string & filename);
    int setGraph(const string & graphfilename, const string & cutfilename);

    int calculateEdgeLength();
    int calculateEdgeWeight();

    double distance(CTarget * x, CTarget * y);
    double distance(CTarget * x, const SmartGraph::Edge & e, SmartGraph::Node & nx, SmartGraph::Node & ne);
    double distance(CTarget * x, const SmartGraph::Node & n, SmartGraph::Node & nx);
    double distance(CTarget * x, SmartGraph::Node n);

    double calculateBarycenter(CVertex * v, vector<CVertex*> & nei);

    int initialMap(Cut & cuts, map<int, int> & seeds);
    int harmonicMap();

    int traceAllPants(const map<int, int>& seeds, map<int, vector<CVertex*>>& pantss);

    int tracePants(int id, int seed, vector<CVertex*>& pants);

	int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants);	
	int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants, SmartGraph::Edge & e0, SmartGraph::Edge & e1, SmartGraph::Edge & e2);

    int findNeighbors( vector<int> & cut, vector<CVertex*> & vs1, vector<CVertex*> & vs2);

    void test();

    int writeMap(string filename);

private:
    CMesh * mesh;
    CGraph * graph;
	Cut cuts;
};
