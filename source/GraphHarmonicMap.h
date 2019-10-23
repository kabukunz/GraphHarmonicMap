#pragma once

#include "Mesh/mesh.h"
#include "Mesh/dynamicmesh.h"
#include "Mesh/boundary.h"
#include "Graph.h"

//using namespace MeshLib;

#ifndef EPS
#define EPS 1e-6
#endif // !EPS

#define ADD_PROPERTY(T, x) \
private:\
    T m_##x; \
public:\
    inline T & x() { return m_##x; } \

class CTarget;

class CHVertex
{
    ADD_PROPERTY(int, index)
    ADD_PROPERTY(CTarget*, target)
    ADD_PROPERTY(bool, fixed)
    ADD_PROPERTY(bool, critical)
    ADD_PROPERTY(bool, critical2)
    ADD_PROPERTY(double, x)
    ADD_PROPERTY(double, y)
    ADD_PROPERTY(bool, cut)
    ADD_PROPERTY(int, cut_id)
    ADD_PROPERTY(bool, cut2)
    ADD_PROPERTY(int, pants)
    ADD_PROPERTY(double, u)
    //ADD_PROPERTY(bool, boundary)

    ADD_PROPERTY(int, nn)
    ADD_PROPERTY(CTarget**, neit)
    ADD_PROPERTY(double*, ew)
    ADD_PROPERTY(double, ewsum)
    ADD_PROPERTY(double, lambda)
    ADD_PROPERTY(double*, bx)
    ADD_PROPERTY(short*, be)
    ADD_PROPERTY(double*, bp)
    ADD_PROPERTY(CTarget**, vx)
    ADD_PROPERTY(double*, fm)
};
class CHEdge
{
    ADD_PROPERTY(double, weight)
    ADD_PROPERTY(int, index)
};
class CHFace
{
    ADD_PROPERTY(int, index)
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



typedef map<int, vector<CVertex*>> Cut;
typedef map<int, set<int>> CutMaps;
typedef map<int, double> WeightMaps;
typedef map<int, CVertex*> Seed;
typedef map<int, vector<CVertex*>> Pants;

class CGraphHarmonicMap
{
public:

public:
    CGraphHarmonicMap();
    ~CGraphHarmonicMap();

    int setMesh(const string & filename);
    int setGraph(const string & graphfilename);

    int calculateEdgeLength();
    int calculateEdgeWeight();

    int runRicciFlow();

    double distance(CTarget * x, CTarget * y);
    double distance(CTarget * x, const SmartGraph::Edge & e, SmartGraph::Node & nx, SmartGraph::Node & ne);
    double distance(CTarget * x, const SmartGraph::Node & n, SmartGraph::Node & nx);
    double distance(CTarget * x, SmartGraph::Node n);

    double calculateBarycenter(CVertex * v);
    double calculateHarmonicEnergy(CVertex * v);
    double calculateHarmonicEnergy();

    int initialMap(string method = string("init"));
    int harmonicMap();

    int traceAllPants();
    int tracePants(int id, CVertex * seed, vector<CVertex*> & pants);

    int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants);
    int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants, SmartGraph::Edge & e0, SmartGraph::Edge & e1, SmartGraph::Edge & e2);
    int embedPants(SmartGraph::Node & node, vector<CVertex*> & pants, vector<SmartGraph::Edge> & edges);

    int findNeighbors(vector<CVertex*> & cut, vector<CVertex*> & vs1, vector<CVertex*> & vs2);
    //int findMeshBoundaries();

    void test();

    int traceCriticalTrajectory();

    int decompose();

    bool hasCriticalPoint(CVertex * v1, CVertex * v2);
    bool hasCriticalPoint(CEdge * e);
    bool hasCriticalPoint(CFace * f);
    CVertex * locateCriticalPoint(CFace * f);
    CVertex * locateCriticalPoint(CEdge * e);

    int output(string filename);
    int outputGraph(const string & filename);

private:
    CMesh * mesh;
    CGraph * graph;
    Cut cuts;
    Seed seeds;
    Pants pantss;
    CutMaps cms;
    WeightMaps wms;

    CDynamicMesh * dmesh;
    CBoundary * boundary;
};
