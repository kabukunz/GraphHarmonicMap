#include "Graph.h"
#include <lemon/dijkstra.h>

CGraph::CGraph() : nodePosition(g), edgeLength(g), dist(g)
{
}


CGraph::~CGraph()
{
}

void CGraph::read(string filename)
{
    g.clear();
    std::ifstream infile(filename);
    string line;

    while (getline(infile, line))
    {
        std::istringstream iss(line);
        string type;
        int id;
        double x, y;
        int i, j;
        iss >> type;
        if (type[0] == '#') continue;

        if (type == "Node")
        {
            iss >> id >> x >> y;
            auto n = g.addNode();
            nodePosition[n] = CPoint2(x, y);
            nodeMap[id] = n;
        }
        if (type == "Arc")
        {
            double l;
            iss >> id >> i >> j >> l;
            auto arc = g.addEdge(nodeMap[i], nodeMap[j]);
            edgeLength[arc] = l;
        }
    }
}

void CGraph::write(string filename)
{
    std::ofstream outfile(filename);
    for (SmartGraph::NodeIt n(g); n != INVALID; ++n)
    {
        outfile << g.id(n) << ": ";
        for (SmartGraph::OutArcIt a(g, n); a != INVALID; ++a)
        {
            outfile << g.id(g.target(a)) << ", ";
        }
        outfile << endl;
    }
}

double CGraph::distance(const SmartGraph::Node & n1, const SmartGraph::Node & n2)
{
    if (nodeDistance.size() == 0) calculateNodeDistance();
    return nodeDistance[make_pair(n1, n2)];
}

void CGraph::calculateNodeDistance()
{
    auto dij = dijkstra(g, edgeLength).distMap(dist);
    for (SmartGraph::NodeIt n1(g); n1 != INVALID; ++n1)
    {
        for (SmartGraph::NodeIt n2(g); n2 != INVALID; ++n2)
        {
            double d = 0.0;
            if (n1 != n2) dij.dist(d).run(n1, n2);
            nodeDistance[make_pair(n1, n2)] = d;
        }
    }
}


