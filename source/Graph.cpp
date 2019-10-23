#include "Graph.h"
#include <lemon/dijkstra.h>

CGraph::CGraph() : edgeLength(g), edgeSign(g), nodeValence(g)
{
}


CGraph::~CGraph()
{
    for (int i = 0; i < g.nodeNum(); ++i)
    {
        delete[] nodeDist[i];
    }
    delete[] nodeDist;
}


SmartGraph::Edge CGraph::addEdge(int i, int j, double length)
{
    while (g.nodeNum() < std::max(i, j) + 1)
    {
        SmartGraph::Node node = g.addNode();
        nodeValence[node] = 0;
    }
    SmartGraph::Node ni = g.nodeFromId(i);
    nodeMap[i] = ni;

    SmartGraph::Node nj = g.nodeFromId(j);
    nodeMap[j] = nj;

    SmartGraph::Edge edge = g.addEdge(ni, nj);
    edgeLength[edge] = length;
    nodeValence[ni] += 1;
    nodeValence[nj] += 1;

    return edge;
}

ostream& operator<<(ostream& os, CGraph& graph)
{
    auto& g = graph.g;
    for (SmartGraph::NodeIt n(g); n != INVALID; ++n)
    {
        os << g.id(n) << ": ";
        for (SmartGraph::OutArcIt a(g, n); a != INVALID; ++a)
        {
            os << g.id(g.target(a)) << ", ";
        }
        os << endl;
    }
    return os;
}

double CGraph::distance(const SmartGraph::Node & n1, const SmartGraph::Node & n2)
{
    int i = g.id(n1);
    int j = g.id(n2);
    return nodeDist[i][j];
}

void CGraph::calculateNodeDistance()
{
    SmartGraph::NodeMap<double> dist(g);
    auto dij = dijkstra(g, edgeLength).distMap(dist);
    nodeDist = new double*[g.nodeNum()];
    for (int i = 0; i < g.nodeNum(); ++i)
    {
        nodeDist[i] = new double[g.nodeNum()];
    }
    for (SmartGraph::NodeIt n1(g); n1 != INVALID; ++n1)
    {
        int i = g.id(n1);
        for (SmartGraph::NodeIt n2(g); n2 != INVALID; ++n2)
        {
            double d = 0.0;
            int j = g.id(n2);
            if (n1 != n2) dij.dist(d).run(n1, n2);
            nodeDist[i][j] = d;
        }
    }
}

int CGraph::colorize()
{
    for (SmartGraph::NodeIt n(g); n != INVALID; ++n)
    {
        for (SmartGraph::OutArcIt oa(g, n); oa != INVALID; ++oa)
        {
            SmartGraph::Edge e = oa;
            edgeSign[e] = 0;
        }
    }
    int s = 1;
    for (SmartGraph::NodeIt n(g); n != INVALID; ++n)
    {
        for (SmartGraph::OutArcIt oa(g, n); oa != INVALID; ++oa)
        {
            SmartGraph::Edge e = oa;
            int sign = edgeSign[e];
            if (sign != 0) s = -sign;
        }
        for (SmartGraph::OutArcIt oa(g, n); oa != INVALID; ++oa)
        {
            SmartGraph::Edge e = oa;
            int sign = edgeSign[e];
            if (sign == 0)
            {
                edgeSign[e] = s;
                s = -s;
            }
        }
    }
    return 0;
}

int CGraph::saveGraph(const string & filename, const map<int, vector<int>> & cuts, const map<int, vector<int>> & pants)
{
    return 0;
}
