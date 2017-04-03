#pragma once

#include "Mesh/mesh.h"
#include "Mesh/dynamicmesh.h"
#include "Mesh/boundary.h"
#include "Graph.h"

//using namespace MeshLib;

#ifndef EPS
#define EPS 1e-7
#endif // !EPS

#define ADD_PROPERTY(T, x) \
private:\
    T m_##x; \
public:\
    inline T & x() { return m_##x; } \

class CTarget;

class CHVertex
{
    ADD_PROPERTY(CTarget*, target)
    ADD_PROPERTY(bool, fixed)
    ADD_PROPERTY(bool, critical)
    ADD_PROPERTY(bool, critical2)
    ADD_PROPERTY(double, x)
    ADD_PROPERTY(double, y)
    ADD_PROPERTY(bool, cut)
    ADD_PROPERTY(bool, cut2)
    ADD_PROPERTY(int, pants)
    ADD_PROPERTY(double*, bx)
};
class CHEdge
{
    ADD_PROPERTY(double, weight)
};
class CHFace
{

};
class CHHalfEdge
{

};
using CMesh = MeshLib::CBaseMesh<CHVertex, CHEdge, CHFace, CHHalfEdge>;
using CDynamicMesh = MeshLib::CDynamicMesh<CHVertex, CHEdge, CHFace, CHHalfEdge>;
using CBoundary = MeshLib::CBoundary<CHVertex, CHEdge, CHFace, CHHalfEdge>;
using CLoop = MeshLib::CLoop<CHVertex, CHEdge, CHFace, CHHalfEdge>;
using CPoint = MeshLib::CPoint;

using CVertex = CMesh::CVertex;
using CEdge = CMesh::CEdge;
using CFace = CMesh::CFace;
using CHalfEdge = CMesh::CHalfEdge;

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
    
public:
    CGraphHarmonicMap();
    ~CGraphHarmonicMap();

    int setMesh(const string & filename);
    int setGraph(const string & graphfilename, const string & cutfilename);

    int calculateEdgeLength();
    int calculateEdgeWeight();

    int runRicciFlow();

    double distance(CTarget * x, CTarget * y);
    double distance(CTarget * x, const SmartGraph::Edge & e, SmartGraph::Node & nx, SmartGraph::Node & ne);
    double distance(CTarget * x, const SmartGraph::Node & n, SmartGraph::Node & nx);
    double distance(CTarget * x, SmartGraph::Node n);

    double calculateBarycenter(CVertex * v, vector<CHalfEdge*> & nei);

    int initialMap(string method = string("init"));
    int harmonicMap();

    int traceAllPants();
    int tracePants(int id, int seed, vector<CVertex*>& pants);

    int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants);
    int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants, SmartGraph::Edge & e0, SmartGraph::Edge & e1, SmartGraph::Edge & e2);

    int findNeighbors(vector<int> & cut, vector<CVertex*> & vs1, vector<CVertex*> & vs2);

    void test();

    int colorizeGraph();

    int traceCriticalTrajectory();

    int decompose();

    bool hasCriticalPoint(CVertex * v1, CVertex * v2);
    bool hasCriticalPoint(CEdge * e);
    bool hasCriticalPoint(CFace * f);
    CVertex * locateCriticalPoint(CFace * f);
    CVertex * locateCriticalPoint(CEdge * e);

    int output(string filename);

private:
    CMesh * mesh;
    CGraph * graph;
    Cut cuts;
    Seed seeds;
    Pants pantss;

    CDynamicMesh * dmesh;
    CBoundary * boundary;
};
