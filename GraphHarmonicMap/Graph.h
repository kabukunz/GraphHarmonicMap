#pragma once

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include "lemon/smart_graph.h"
#include "lemon/static_graph.h"
#include "Geometry/Point2.h"

using namespace std;
using namespace lemon;
using namespace MeshLib;

class CGraph
{
public:
    CGraph();
    ~CGraph();

    void read(string filename);
    void write(string filename);

	double distance(const SmartGraph::Node & n1, const SmartGraph::Node & n2);

    SmartGraph g;
    SmartGraph::NodeMap<CPoint2> nodePosition;
    SmartGraph::EdgeMap<double> edgeLength;
	SmartGraph::NodeMap<double> dist;

private:
	map<int, SmartGraph::Node> nodeMap;
	map<pair<SmartGraph::Node, SmartGraph::Node>, double> nodeDistance;
	void calculateNodeDistance();
};
