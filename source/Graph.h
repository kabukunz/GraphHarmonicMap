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

    int read(string filename);
    int write(string filename);

    void calculateNodeDistance();
    double distance(const SmartGraph::Node & n1, const SmartGraph::Node & n2);

    SmartGraph g;
    SmartGraph::EdgeMap<double> edgeLength;
    SmartGraph::EdgeMap<int> edgeSign;
    SmartGraph::NodeMap<int> nodeValence;

private:
    map<int, SmartGraph::Node> nodeMap;
    double ** nodeDist;
    
};
