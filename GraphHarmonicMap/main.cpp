// GraphHarmonicMap.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GraphHarmonicMap.h"

int main(int argc, char * argv[])
{
    if (argc < 3)
    {
        cout << "usage: GraphHarmonicMap meshfilename graphfilename [cutfilename]" << endl;
        return -1;
    }
    CGraphHarmonicMap * map = new CGraphHarmonicMap();
    string meshfilename(argv[1]);
    string graphfilename(argv[2]);
    string cutfilename(argv[3]);

    map->setMesh(meshfilename);
    map->setGraph(graphfilename, cutfilename);
    map->harmonicMap();
    //map->test();

    map->writeMap(meshfilename);

    return 0;
}

