// GraphHarmonicMap.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GraphHarmoicMap.h"

int main(int argc, char * argv[])
{
    if (argc < 3)
    {
        cout << "usage: GraphHarmonicMap meshfilename graphfilename [cutfilename]" << endl;
        return -1;
    }
    CGraphHarmoicMap * map = new CGraphHarmoicMap();
    string meshfilename(argv[1]);
    string graphfilename(argv[2]);
    string cutfilename(argv[3]);

    map->setMesh(meshfilename);
    map->setGraph(graphfilename, cutfilename);
    map->harmonicMap();

    map->writeMap(meshfilename);

    return 0;
}

