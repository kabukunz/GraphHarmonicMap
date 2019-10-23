#pragma once

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include "lemon/smart_graph.h"
#include "lemon/static_graph.h"

using namespace std;
using namespace lemon;

class CGraph
{
public:
    CGraph();
    ~CGraph();

    SmartGraph::Edge addEdge(int i, int j, double length);

    void calculateNodeDistance();
    double distance(const SmartGraph::Node & n1, const SmartGraph::Node & n2);

    int colorize();
    int saveGraph(const string & filename, const map<int, vector<int>> & cuts, const map<int, vector<int>> & pants);

    SmartGraph g;
    SmartGraph::EdgeMap<double> edgeLength;
    SmartGraph::EdgeMap<int> edgeSign;
    SmartGraph::NodeMap<int> nodeValence;

    friend ostream& operator<<(ostream& os, CGraph& graph);

private:
    map<int, SmartGraph::Node> nodeMap;
    double ** nodeDist;

};
