#pragma once

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include "lemon/smart_graph.h"
#include "lemon/static_graph.h"
#include "Mesh/Base.h"

using namespace std;
using namespace lemon;

class CGraph : public CBase
{
public:
    CGraph();
    ~CGraph();

    void read(string filename);
    void write(string filename);

    double distance(const SmartGraph::Node & n1, const SmartGraph::Node & n2);

	SmartGraph g;
	SmartGraph::EdgeMap<double> edgeLength;
	SmartGraph::EdgeMap<int> edgeSign;
	SmartGraph::NodeMap<double> dist;

private:
    map<int, SmartGraph::Node> nodeMap;
    map<pair<SmartGraph::Node, SmartGraph::Node>, double> nodeDistance;
    void calculateNodeDistance();
};
