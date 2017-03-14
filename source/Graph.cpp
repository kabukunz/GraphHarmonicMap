#include "Graph.h"
#include <lemon/dijkstra.h>

CGraph::CGraph() : edgeLength(g), edgeSign(g), nodeValence(g)
{
}


CGraph::~CGraph()
{
}

int CGraph::read(string filename)
{
    g.clear();
    std::ifstream infile(filename);
    if (!infile.good())
    {
        cout << "can't open graph file: " << filename << endl;
        exit(-1);
        return -1;
    }
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
            nodeMap[id] = n;
            nodeValence[n] = 0;
        }
        if (type == "Arc")
        {
            double l;
            iss >> id >> i >> j >> l;
            auto arc = g.addEdge(nodeMap[i], nodeMap[j]);
            edgeLength[arc] = l;
            nodeValence[nodeMap[i]] += 1;
            nodeValence[nodeMap[j]] += 1;
        }
    }
    return 0;
}

int CGraph::write(string filename)
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
    return 0;
}

double CGraph::distance(const SmartGraph::Node & n1, const SmartGraph::Node & n2)
{
    if (nodeDistance.size() == 0) calculateNodeDistance();
    return nodeDistance[make_pair(n1, n2)];
}

void CGraph::calculateNodeDistance()
{
    SmartGraph::NodeMap<double> dist(g);
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


