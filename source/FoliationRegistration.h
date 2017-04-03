#pragma once

#include "Mesh/mesh.h"
#include "Mesh/iterators.h"
#include "Graph.h"

using namespace MeshLib;

class CFoliationRegistration
{
public:
    CFoliationRegistration();
    ~CFoliationRegistration();

    int setSource(CMesh * mesh, CGraph * graph);
    int setTarget(CMesh * mesh, CGraph * graph);
    int initialMap();
    int diffuse();

    int output();

private:
    CMesh * sourceMesh;
    CMesh * targetMesh;
    CGraph * sourceGraph;
    CGraph * targetGraph;
};

