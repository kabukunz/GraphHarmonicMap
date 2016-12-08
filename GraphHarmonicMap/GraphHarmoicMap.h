#pragma once

#include "Mesh/Mesh.h"
#include "Mesh/Iterators.h"
#include "Graph.h"

using namespace MeshLib;

/*
* represent a target on the graph, identified by an arc, and the direction, and the length from starting point
* for a loop, there is only one starting point, but two directions: +1, -1;
* for otherwise, there are two starting points +1 and -1, direction +1 means starting from +1, -1 starting from -1
*/
class CTarget
{
public:
    SmartGraph::Edge edge;
    int direction;
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

    int calculateBarycenter(CVertex* v);
    int harmonicMap();

    int writeMap(string filename);

private:
    CMesh * mesh;
    CGraph * graph;
};
