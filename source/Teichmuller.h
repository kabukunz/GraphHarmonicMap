#pragma once
#include "Mesh/Mesh.h"
#include "Mesh/Iterators.h"
#include "Graph.h"

using namespace MeshLib;

class CTeichmuller
{
public:
	CTeichmuller();
	~CTeichmuller();

private:
	CMesh * mesh;
};

