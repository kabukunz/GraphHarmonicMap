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


    SmartGraph g;
    SmartGraph::NodeMap<CPoint2> nodePosition;
    SmartGraph::EdgeMap<double> edgeLength;

private:
	map<int, SmartGraph::Node> nodeMap;
};
