#pragma once

#include "Mesh/Mesh.h"
#include "Mesh/Iterators.h"
#include "Graph.h"

using namespace MeshLib;

#ifndef EPS
#define EPS 1e-6
#endif // !EPS


/*
* represent a target on the graph, identified by an edge, and the starting node, and the length from starting node
*/
class CTarget
{
public:
    SmartGraph::Edge edge;
    SmartGraph::Node node;
    double length;
};

typedef map<int, vector<int>> Cut;
typedef map<int, int> Seed;
typedef map<int, vector<CVertex*>> Pants;

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

    int initialMap();
    int harmonicMap();

    int traceAllPants();

    int tracePants(int id, int seed, vector<CVertex*>& pants);

    int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants);
    int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants, SmartGraph::Edge & e0, SmartGraph::Edge & e1, SmartGraph::Edge & e2);

    int findNeighbors( vector<int> & cut, vector<CVertex*> & vs1, vector<CVertex*> & vs2);

    int colorizeGraph();
    int traceCriticalTrajectory();

    void test();

    int writeMap(string filename);

private:
    CMesh * mesh;
    CGraph * graph;
    Cut cuts;
    Seed seeds;
    Pants pantss;
};
