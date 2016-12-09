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

class CGraphHarmoicMap
{
public:
    CGraphHarmoicMap();
    ~CGraphHarmoicMap();

    int setMesh(string filename);
    int setGraph(string filename);

    int calculateEdgeLength();
    int calculateEdgeWeight();

	double distance(CTarget * x, CTarget * y);
	double distance(CTarget * x, const SmartGraph::Edge & e, SmartGraph::Node & nx, SmartGraph::Node & ne);	
	double distance(CTarget * x, const SmartGraph::Node & n, SmartGraph::Node & nx);
	double distance(CTarget * x, SmartGraph::Node n);

    double calculateBarycenter(CVertex* v);
    int harmonicMap();

	void test();

    int writeMap(string filename);

private:
    CMesh * mesh;
    CGraph * graph;
};
