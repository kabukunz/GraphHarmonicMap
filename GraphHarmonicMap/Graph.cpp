#include "Graph.h"


CGraph::CGraph() : nodePosition(g), edgeLength(g)
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


